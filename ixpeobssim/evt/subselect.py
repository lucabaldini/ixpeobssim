#!/usr/bin/env python
#
# Copyright (C) 2015--2022, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""Event list filtering facilities.
"""

from __future__ import print_function, division

import numpy
import regions

from ixpeobssim.evt.event import xEventFile
from ixpeobssim.instrument.gpd import GPD_DEFAULT_FIDUCIAL_HALF_SIDE_X,\
    GPD_DEFAULT_FIDUCIAL_HALF_SIDE_Y
from ixpeobssim.instrument.mma import fiducial_backscal
from ixpeobssim.utils.astro import angular_separation, ds9_region_filter_sky
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.profile import timing
from ixpeobssim.utils.units_ import degrees_to_arcmin, arcmin_to_arcsec, arcmin_to_degrees


# pylint: disable=invalid-name, too-many-locals, consider-using-f-string


class xEventSelect:

    """Base class for event subselection.
    """

    def __init__(self, file_path, **kwargs):
        """Constructor.
        """
        self.event_file = xEventFile(file_path)
        self.kwargs = kwargs
        self._process_kwargs()

    def _process_kwargs(self):
        """Check the keyword arguments.
        """
        if self.get('suffix') is None:
            suffix = 'select'
        else:
            suffix = self.get('suffix')
        evfile = self.event_file.file_path()
        outfile = evfile.replace('.fits', '_%s.fits' % suffix)
        self.set('outfile', outfile)
        if self._acceptance_cone_required():
            obj_ra, obj_dec = self.event_file.wcs_reference()
            if self.get('ra') is None:
                if obj_ra is None:
                    abort('Please set the --ra command-line switch explicitely')
                self.set('ra', obj_ra)
            if self.get('dec') is None:
                if obj_dec is None:
                    abort('Please set the --dec command-line switch explicitely')
                self.set('dec', obj_dec)

    def get(self, key, default=None):
        """Convenience method to address the keyword aguments.
        """
        return self.kwargs.get(key, default)

    def set(self, key, value):
        """Convenience method to set keyword arguments.
        """
        logger.info('Setting %s to %s...', key, value)
        self.kwargs[key] = value

    def _acceptance_cone_required(self):
        """Convenience method to find out whether acceptance cone is actually
        required.

        This is true if we are using the --rad or --innerrad command-line switches.
        """
        return (self.get('rad'), self.get('innerrad')) !=  (None, None)

    def time_selected(self):
        """Return True if selecting in time, i.e., tmin and/or tmax are not None.
        """
        return (self.get('tmin'), self.get('tmax')) != (None, None)

    def phase_selected(self):
        """Return True if selecting in phase, i.e., phasemin and/or phasemax are not None.
        """
        return (self.get('phasemin'), self.get('phasemax')) != (None, None)

    def mask_selected(self):
        """Return True if selecting with event mask directly,
        """
        return self.get('mask') is not None

    def _validate(self):
        """Make sure that the event selection is valid.

        We want to make sure that any selection in time and phase is within the
        physical bounds to make it harder to make mistakes when updating the
        header keywords in the output files.

        Also, we are not supporting simultaneous selections in time and phase.

        Note this is aborting if any of the relevant conditions is not met.
        """
        if self.time_selected() and self.phase_selected():
            abort('Combined selections on time and phase not currently supported')
        if self._acceptance_cone_required() and self.get('regfile') is not None:
            abort('Combined circle/annulus and ds9 regions not currently supported')
        # Check the time-related keywords.
        primary_header = self.event_file.primary_header
        tstart, tstop = primary_header.get('TSTART'), primary_header.get('TSTOP')
        # Make sure that TSTART <= tmin <= TSTOP...
        tmin = self.get('tmin')
        if tmin is not None and (tmin < tstart or tmin > tstop):
            abort('tmin (%.3f) outside the physical bounds [%.3f--%.3f]' % (tmin, tstart, tstop))
        # Make sure that TSTART <= tmax <= TSTOP...
        tmax = self.get('tmax')
        if tmax is not None and (tmax < tstart or tmax > tstop):
            abort('tmax (%.3f) outside the physical bounds [%.3f--%.3f]' % (tmax, tstart, tstop))
        # Make sure that tmin <= tmax
        if tmin is not None and tmax is not None and tmin >= tmax:
            abort('tmin (%.3f) >= tmax (%.3f)' % (tmin, tmax))
        # Make sure 0. <= phasemin <= 1.
        phasemin = self.get('phasemin')
        if phasemin is not None and (phasemin < 0. or phasemin > 1.):
            abort('phasemin (%.3f) outside the physical bounds [0.--1.]' % phasemin)
        # Make sure 0. <= phasemax <= 1.
        phasemax = self.get('phasemax')
        if phasemax is not None and (phasemax < 0. or phasemax > 1.):
            abort('phasemax (%.3f) outside the physical bounds [0.--1.]' % phasemax)
        # Make sure that phasemin <= phasemax
        if phasemin is not None and phasemax is not None and phasemin >= phasemax:
            abort('phasemin (%.3f) >= phasemax (%.3f)' % (phasemin, phasemax))

    def _time_header_keywords(self, mask, livetime_ratio_threshold=0.02):
        """Calculate the updated time header keywords to be written in the output
        file downstream the selection.

        This is actually one of the most difficult thing to get right, and it is
        paramount for the filtered files to operate with all the analysis ecosystem
        in just the same way the original files were doing.

        .. warning::

           Note that I am afraid multiple selections in phase will be essentially
           impossible to get right.
        """
        # If we are not correcting the livetime, likewise.
        if not self.get('ltimeupdate'):
            return None

        # And now we are in business...
        header_keywords = {}
        primary_header = self.event_file.primary_header
        average_dtpe = self.event_file.average_deadtime_per_event()
        scaled_livetime = None
        # Retrieve the livetime and convert it to s.
        livetime = 1.e-6 * self.event_file.livetime_data()
        total_livetime_sum = livetime.sum()
        logger.info('Initial livetime sum: %.3f s', total_livetime_sum)
        logger.info('Initial number of events: %d s', self.event_file.num_events())
        num_events = mask.sum()
        livetime_sum = livetime[mask].sum()
        logger.info('Residual livetime sum after time/phase selection: %.3f s', livetime_sum)
        logger.info('Remanining events after time/phase selection: %d', num_events)

        # Are we selecting in time? In this case we definitely want to update
        # TSTART, TSTOP and ONTIME, as well as the livetime-related quantities.
        if self.time_selected():
            tstart = self.get('tmin', primary_header['TSTART'])
            tstop = self.get('tmax', primary_header['TSTOP'])
            ontime = tstop - tstart
            selection_frac = ontime / primary_header.get('ONTIME')
            logger.info('Raw time selection fraction: %.5f', selection_frac)
            scaled_livetime = ontime - average_dtpe * num_events
            logger.debug('TSTART: %.3f -> %.3f', primary_header.get('TSTART'), tstart)
            logger.debug('TSTOP : %.3f -> %.3f', primary_header.get('TSTOP'), tstop)
            logger.debug('ONTIME: %.3f -> %.3f', primary_header.get('ONTIME'), ontime)
            header_keywords['TSTART'] = tstart
            header_keywords['TSTOP'] = tstop
            header_keywords['ONTIME'] = ontime

        # Are we selecting on phase? In this case we leave TSTART and TSTOP
        # untouched, but we do want to update the ONTIME keyword based on the
        # phase selection interval, so that DEADC downstream has a sensible value.
        elif self.phase_selected():
            phasemin, phasemax = self.get('phasemin', 0.), self.get('phasemax', 1.)
            selection_frac = phasemax - phasemin
            logger.info('Raw phase selection fraction: %.5f', selection_frac)
            ontime = primary_header.get('ONTIME') * selection_frac
            scaled_livetime = ontime - average_dtpe * num_events
            logger.debug('ONTIME: %.3f -> %.3f', primary_header.get('ONTIME'), ontime)
            header_keywords['ONTIME'] = ontime

        # This might be helpful for debug. When the scaled livetime and the
        # simple livetime sum are not in agreement with each other, this might
        # indicate that the former is more accurate.
        if scaled_livetime is not None:
            logger.info('Scaled estimated livetime: %.3f', scaled_livetime)
            _ratio = scaled_livetime / livetime_sum
            logger.debug('Livetime estimate ratio: %.5f', _ratio)
            if abs(_ratio - 1) > livetime_ratio_threshold:
                msg = 'The livetime sum and the scaled livetime differ by more than %.2f%%' %\
                    livetime_ratio_threshold
                logger.warning(msg)
                logger.warning('This might indicate that you have to look into it manually...')

        # Finally, the livetime stuff.
        if self.get('ltimealg') == 'LTSUM':
            livetime = livetime_sum
        elif self.get('ltimealg') == 'LTSCALE':
            livetime = scaled_livetime
        header_keywords['LIVETIME'] = livetime
        header_keywords['DEADC'] = livetime / ontime
        logger.debug('Updated header keywords: %s', header_keywords)
        return header_keywords

    def _time_selection_mask(self):
        """Apply the time selection.

        This is done first, and separately from the other selections, as it
        impacts the livetime calculation. The function returns a boolean mask of
        the same length as the input photon list.
        """
        mask = numpy.ones(self.event_file.num_events(), 'bool')
        tmin, tmax = self.get('tmin'), self.get('tmax')
        if tmin is not None:
            mask *= (self.event_file.time_data() >= tmin)
        if tmax is not None:
            mask *= (self.event_file.time_data() < tmax)
        if self.get('tinvert'):
            mask = numpy.logical_not(mask)
        return mask

    def _phase_selection_mask(self):
        """Apply the phase selection.

        This is done first, and separately from the other selections, as it
        impacts the livetime calculation. The function returns a boolean mask of
        the same length as the input photon list.
        """
        mask = numpy.ones(self.event_file.num_events(), 'bool')
        phasemin, phasemax = self.get('phasemin'), self.get('phasemax')
        if phasemin is not None:
            mask *= (self.event_file.phase_data() >= phasemin)
        if phasemax is not None:
            mask *= (self.event_file.phase_data() < phasemax)
        if self.get('phaseinvert'):
            mask = numpy.logical_not(mask)
        return mask

    def _direct_selection_mask(self):
        """Apply the direct selection reading the mask from file
        """
        with open(self.get('mask'), 'rb') as f:
            mask_in = numpy.load(f)
        return mask_in

    @staticmethod
    def _default_fiducial_backscal():
        """Return the default fiducial backscal for the entire detector.
        """
        return fiducial_backscal(GPD_DEFAULT_FIDUCIAL_HALF_SIDE_X, GPD_DEFAULT_FIDUCIAL_HALF_SIDE_Y)

    @staticmethod
    def _annulus_backscal(inner_radius=None, outer_radius=None):
        """Calculate the value to be written in the BACKSCAL header keyword for
        a given annulus selection.

        This is essentially the area of the annulus in arcsec squared.

        Arguments
        ---------
        inner_radius : float
            The inner radius of the annulus in arcminutes.

        outer_radius : float
            The outer radius of the annulus in arcminutes.

        Return
        ------
        The area of the annulus in arcsec squared.
        """
        if inner_radius is None:
            inner_area = 0.
        else:
            inner_area = numpy.pi * arcmin_to_arcsec(inner_radius)**2.
        if outer_radius is None:
            outer_area = xEventSelect._default_fiducial_backscal()
        else:
            outer_area = numpy.pi * arcmin_to_arcsec(outer_radius)**2.
        return outer_area - inner_area

    @timing
    def _region_backscal(self, region_file_path, num_samples=10000000, half_side=10.):
        """Calculate the value to be written in the BACKSCAL header keyword for
        an arbitrary region shape passed via DS9.

        The procedure is a simple implementation of an hit and miss algorithm,
        where we throw a bunch of sky coordinates uniformly distributed within
        a square in the sky and count the fraction of them falling inside the
        region (or the logical or of the regions).

        The 10' hald_side is necessary to contain all the field of view including
        DU rotation and dithering.

        Arguments
        ---------
        region_file_path : str
            The path to the ds9 region file.

        num_samples : int
            The number of samples to be extracted for the hit or miss integration

        half_side : float
            The half-side of the sampling region in arcmin.
        """
        # pylint: disable=protected-access
        # Remember: we measure the area in arcsec^2, and we throw sky coordinates
        # in degrees.
        total_area = 4. * arcmin_to_arcsec(half_side)**2.
        half_side = arcmin_to_degrees(half_side)
        # Throw a bunch of random coordinates in the sky, sampling a square region
        # centered on the reference of the underlying WCS and larger than the
        # full FOV of the instrument.
        ra0, dec0 = self.event_file.wcs_reference()
        delta_ra = numpy.random.uniform(-half_side, half_side, size=num_samples)
        ra = ra0 + delta_ra / numpy.cos(numpy.radians(dec0))
        dec = numpy.random.uniform(dec0 - half_side, dec0 + half_side, size=num_samples)
        # Read the regions. Note that, for the actual event mask downstream,
        # the relevant call is to xEventFile.ds9_region_file_mask(), which takes the
        # logical or of the regions (when the file contains more than one of them)
        # and here we are doing exactly the same thing.
        regs = regions.Regions.read(region_file_path)
        # Count the events inside the region.
        hit = ds9_region_filter_sky(ra, dec, self.event_file._wcs, *regs).sum()
        backscal = total_area * hit / num_samples
        logger.info ('Sampling area: %s, %d points / %d hit', total_area, hit, num_samples)
        logger.info ('BACKSCAL value: %.3f arcsec^2', backscal)
        return total_area * hit / num_samples


    def select(self):
        """Select the events and write the output file.
        """
        # pylint: disable=too-many-branches, too-many-statements
        self._validate()
        # Prepare a list of the header keywords that needs to be updated in the
        # filtered file---these will be passed to xEventFile.write_fits_selected()
        # in order for the keywords to be written (or updated) in the 'Primary',
        # 'EVENTS', and 'GTI' extensions
        header_keywords = {}

        logger.info('Running event selection with kwargs %s...', self.kwargs)
        total_num_events = self.event_file.num_events()
        if self.time_selected():
            mask = self._time_selection_mask()
            _keywords = self._time_header_keywords(mask)
            if _keywords is not None:
                header_keywords.update(_keywords)
        elif self.phase_selected():
            mask = self._phase_selection_mask()
            _keywords = self._time_header_keywords(mask)
            if _keywords is not None:
                header_keywords.update(_keywords)
        elif self.mask_selected():
            mask = self._direct_selection_mask()
        else:
            mask = numpy.ones(total_num_events, dtype=bool)

        # Apply the various other selections one by one: energy...
        emin, emax = self.get('emin'), self.get('emax')
        energy_mask = numpy.ones(total_num_events, dtype=bool)
        if emin is not None:
            energy_mask *= (self.event_file.energy_data(self.get('mc')) >= emin)
        if emax is not None:
            energy_mask *= (self.event_file.energy_data(self.get('mc')) < emax)
        if self.get('einvert'):
            energy_mask = numpy.logical_not(energy_mask)
        mask *= energy_mask

        # ...acceptance cone...
        if self._acceptance_cone_required():
            ra, dec = self.event_file.sky_position_data(self.get('mc'))
            ra0, dec0 = self.get('ra'), self.get('dec')
            sep = degrees_to_arcmin(angular_separation(ra, dec, ra0, dec0))
            rad, innerrad = self.get('rad'), self.get('innerrad')
            if rad is not None:
                mask *= (sep <= rad)
            if innerrad is not None:
                mask *= (sep >= innerrad)
            header_keywords['BACKSCAL'] = self._annulus_backscal(innerrad, rad)

        # ... if a region file is passed, go ahead and process the region list...
        # Note that we do not allow simultaneous selection with --rad, --innerrad
        # and ds9 region files.
        regfile = self.get('regfile')
        if regfile is not None:
            reg_mask = self.event_file.ds9_region_file_mask(regfile, mc=self.get('mc'))
            backscal = self._region_backscal(regfile)
            if self.get('reginvert'):
                reg_mask = numpy.logical_not(reg_mask)
                backscal = xEventSelect._default_fiducial_backscal() - backscal
            mask *= reg_mask
            header_keywords['BACKSCAL'] = backscal

        # ... and source ID.
        for srcid in self.get('mcsrcid'):
            mask *= (self.event_file.srcid_data() == srcid)

        # Ready to write the output file!
        logger.info('Done, %d events out of %d left.', mask.sum(), total_num_events)
        history = 'Processed with xpselect with arguments %s.' % self.kwargs
        self.event_file.write_fits_selected(mask, self.get('outfile'), header_keywords,
            history, self.get('overwrite'))
        return self.get('outfile')

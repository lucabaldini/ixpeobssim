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

"""Module encapsulating the event structure and related facilities.
"""

from __future__ import print_function, division

import numbers

from astropy.io import fits
from astropy import wcs
import numpy

from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.core.fitsio import FITS_TO_NUMPY_TYPE_DICT, set_tlbounds
from ixpeobssim.core.hist import xHistogram2d, xScatterPlot
from ixpeobssim.evt.fmt import set_telescope_header_keywords, set_time_header_keywords,\
    set_object_header_keywords, set_version_keywords
from ixpeobssim.evt.fmt import xLvl2PrimaryHDU, xBinTableHDUEvents, xBinTableHDUMonteCarlo,\
    xBinTableHDUGTI, xBinTableHDURoiTable, xBinTableHDUSpacecraftData, xBinTableHDUTimeline,\
    xBinTableHDUOCTI
from ixpeobssim.evt.fmt import _TIME_HEADER_KEYWORDS, WCS_ORIGIN, _SKYCOORD_NUM_SIDE_PIXELS
from ixpeobssim.evt.fmt import standard_xy_to_radec, build_standard_wcs
from ixpeobssim.evt.kislat2015 import xStokesAnalysis
from ixpeobssim.instrument.charging import xEnergyFluxCube, read_charging_parameters
from ixpeobssim.instrument.charging import read_charging_map, create_charging_map_extension
from ixpeobssim.instrument.gpd import within_fiducial_rectangle
from ixpeobssim.instrument.mma import parse_dithering_kwargs
from ixpeobssim.irf.ebounds import NUM_CHANNELS, channel_to_energy, TLMIN, TLMAX
from ixpeobssim.utils.astro import read_ds9, ds9_region_filter_sky, angular_separation
from ixpeobssim.utils.os_ import check_input_file
from ixpeobssim.utils.profile import timing
from ixpeobssim.utils.time_ import current_datetime_string_utc
from ixpeobssim.utils.time_ import met_to_string_utc
from ixpeobssim.utils.units_ import degrees_to_arcmin


# pylint: disable=invalid-name, too-many-ancestors, too-many-arguments, no-member
# pylint: disable=too-many-public-methods, too-many-locals, too-many-lines, too-many-statements



class xBaseEventList(dict):

    """Base class for event lists.

    This is the base class for objects that is created and filled at simulation time
    and then written to file in FITS format. It is essentially a dictionary
    containing all the column data for the output FITS file, providing
    additional facilities for addding, sorting in time and trimming event lists.

    In order to be robust against name changes, subclasses should provide
    all the necessary interfaces to retrieve and set the column values
    (i.e., if the name of any of the columns changes, in principle it should
    not be necessary to make modifications beyond this file). Note that the
    colums related to times and source ID are encapsulated into this base class.

    If the constructor is called with no argument an empty event list is
    created. Alternatively the user can initilize the event list with an array
    of event times (we did pick the time because typically this is the first
    dynamical variable being calculated). In that case the TRG_ID column is also
    filled right away and, if the optional source_id argument is not None, the
    SRC_ID as well. For completeness, this mechanism was added because we
    noticed that there were many calls to xEventList.fill_trigger_id() and
    xEventList.set_source_id() floating around in different places. (Note that,
    altough we need to overwrite the trigger identifier values at the end of the
    loop over the model components, after we have applied the deadtime, we do
    need to fill the corresponding column at the level of the single component
    because all the code downstream, e.g., the xEventList.trim() method, need
    *all* the columns to be present, and have the same length, to operate.)
    """

    _TABLE_CLASSES = ()

    def __init__(self, time_=None, source_id=None):
        """Constructor.
        """
        dict.__init__(self)
        self.specs = sum([cls.spec_names_and_types() for cls in self._TABLE_CLASSES], [])
        self.spec_names = sum([cls.spec_names() for cls in self._TABLE_CLASSES], [])
        for name, dtype in self.specs:
            dtype = FITS_TO_NUMPY_TYPE_DICT[dtype]
            self[name] = numpy.array([], dtype)
        # The class member for the number of events is "private" as all the
        # bookkeeping is handled internally, and we don't want users to change
        # that from outiside the class.
        self._num_events = 0
        # If we are passing the event times to the constructor we do fill the
        # appropriate column, here---and the trigger identifier, too.
        if time_ is not None:
            self.set_time_columns(time_)
        # And if we are passing the source identifier, and assuming the list is
        # not empty, fill that columns too.
        if source_id is not None and self._num_events > 0:
            self.set_source_id(source_id)

    def num_events(self):
        """Return the number of events in the event list.
        """
        return self._num_events

    def empty(self):
        """Return whether the event list is empty.
        """
        return self._num_events == 0

    def _set_column(self, name, data):
        """Set a column array.

        Arguments
        ---------
        name : string
            The name of the column

        data : array or number
            The actual data to put in the column. (If `data` is a number,
            an array of the proper length is automatically created, assuming
            that the `lenght` class member is defined.)
        """
        # Make sure the column name exists.
        assert name in self
        # If data is a number, we want to fill a numpy array with the right
        # type and identical values. (Of course this requires that the event
        # list is not empty, otherwise we have no way of figuring out the length
        # of the array.)
        if isinstance(data, numbers.Number):
            assert self._num_events > 0
            data = numpy.full(self._num_events, data, self[name].dtype)
        # If the event list is empty we know the target length...
        if self.empty():
            self._num_events = len(data)
        # ... otherwise make sure there is no length mismatch.
        elif len(data) != self._num_events:
            logger.error('Addind column %s with %d events to an event list with %d events...',\
                name, len(data), self._num_events)
            abort('Column length mismatch')
        # And now we're good to go---set the actual column data.
        self[name] = data

    def time(self):
        """Return the event times.
        """
        return self.get('TIME')

    @staticmethod
    def split_event_time(time_):
        """Split the event time into its interger and fractional part.

        While in hardware we will be measuring the seconds and microseconds
        separately, and then build the floating-point representation of the
        event time by summing the two, in the simulation realm we do the
        opposite, i.e., we extract the event times as floating-point numbers
        and then reverse-engineer the integral and fractional parts.

        Warning
        -------
        This is a crude implementation of the concept, where we essentially
        take the floor of the integral and fractional parts of the event time
        (modulo the conversion from seconds to microseconds for the second).
        In the future we might make this more complicated if we need to
        simulate the instrument timing in a more detailed fashion.
        """
        seconds = numpy.floor(time_).astype(numpy.int32)
        microseconds = numpy.floor((time_ - seconds) * 1.e6).astype(numpy.int32)
        return (seconds, microseconds)

    @staticmethod
    def array_is_sorted(array):
        """Return whether a given array is sorted.

        This is an interesting topic per se, see
        https://stackoverflow.com/questions/47004506
        I dropped this here assuming that we're not concerned by adding a step
        which is O(n) from an algorithmic standpoint.
        """
        return (numpy.diff(array) >= 0.).all()

    def set_time_columns(self, time_, check_sorted=True):
        """Populate the time-related columns.
        """
        if check_sorted and not self.array_is_sorted(time_):
            abort('Input time array to event list is not sorted %s' % time_)
        seconds, microseconds = self.split_event_time(time_)
        self._set_column('TIME', time_)
        self._set_column('SEC', seconds)
        self._set_column('MICROSEC', microseconds)

    def src_id(self):
        """Return the source ID.
        """
        return self.get('SRC_ID')

    def set_source_id(self, src_id):
        """Set the source ID.
        """
        self._set_column('SRC_ID', src_id)

    def __add__(self, other):
        """Concatenate two event lists.
        """
        # pylint: disable=protected-access
        list_ = self.__class__()
        for name in self.spec_names:
            list_._set_column(name, numpy.append(self[name], other[name]))
        return list_

    def sort(self):
        """Sort the event list based on the event time.
        """
        _index = numpy.argsort(self.time())
        for name in self.spec_names:
            self._set_column(name, self[name][_index])

    def trim(self, mask):
        """Trim all the data columns according to a generic input mask.
        """
        assert len(mask) == self._num_events
        for name in self.spec_names:
            try:
                self[name] = self[name][mask]
            except IndexError as e:
                abort('Cannot trim column %s (%s)' % (name, e))
        self._num_events = numpy.count_nonzero(mask)

    @timing
    def apply_vignetting_base(self, ra, dec, energy, vign, ra_pnt, dec_pnt):
        """Base function to apply the effective area vignetting to the list.

        Note that sub-classes are responsible for providing their own definition
        of the sky coordinates and energy to apply the function.

        .. note::

            Note that angular_separation returns values in degrees, while the
            vignetting function is expressed in terms of the off-axis angle
            expressed in arcminutes.
            See also https://bitbucket.org/ixpesw/ixpeobssim/issues/423
        """
        logger.info('Applying vignetting to the event list...')
        num_events = self.num_events()
        separation = angular_separation(ra, dec, ra_pnt, dec_pnt)
        separation = degrees_to_arcmin(separation)
        vignetting = vign(energy, separation)
        self.trim(numpy.random.random(vignetting.shape) <= vignetting)
        frac = 100. * self.num_events() / float(num_events)
        logger.info('Done, %d/%d (%.2f%%) events remaining.', self.num_events(), num_events, frac)

    @timing
    def apply_fiducial_area_base(self, detx, dety):
        """Trim out the events falling outside the detector active area.

        Note that sub-classes are responsible for providing their own definition
        of the detector coordinates to apply the function.
        """
        logger.info('Applying GPD fiducial area cut to the event list...')
        num_events = self.num_events()
        mask = within_fiducial_rectangle(detx, dety)
        self.trim(mask)
        frac = 100. * self.num_events() / float(num_events)
        logger.info('Done, %d/%d (%.2f%%) events remaining.', self.num_events(), num_events, frac)



class xEventList(xBaseEventList):

    """Class describing an ixpeobssim event list.
    """

    _TABLE_CLASSES = (xBinTableHDUEvents, xBinTableHDUMonteCarlo)

    def __init__(self, time_=None, source_id=None):
        """Constructor.
        """
        xBaseEventList.__init__(self, time_, source_id)
        if time_ is not None:
            self.fill_trigger_id_unphysical()
            self.fill_livetime_unphysical()
            self.set_num_clusters(1)
            self.set_weights(1.)
            self.set_mc_gain_column()
        # We store the final charging status in this list (it's not a single
        # item as we may have two different charging processes going on)
        self.__charging_maps = []

    def livetime(self):
        """Return the event livetimes.
        """
        return self.get('LIVETIME')

    def pi(self):
        """Return the pulse invariant for the events.
        """
        return self.get('PI')

    def mc_energy(self):
        """Return the Monte Carlo event energies.
        """
        return self.get('MC_ENERGY')

    def energy(self):
        """Return the reconstructed event energies.
        """
        return self.get('ENERGY')

    def mc_sky_coordinates(self):
        """Return the sky positions of the events.
        """
        return self.get('MC_RA'), self.get('MC_DEC')

    def phi(self):
        """Return the reconstructed azimuthal angle for the events.
        """
        return self.get('PHI')

    def q(self):
        """Return the Q Stokes parameter.
        """
        return self.get('Q')

    def u(self):
        """Return the U Stokes parameter.
        """
        return self.get('U')

    def stokes_parameters(self):
        """Return the Stokes parameters.
        """
        return self.q(), self.u()

    def detector_coordinates(self):
        """Return the event position in detector coordinates.
        """
        return self.get('DETX'), self.get('DETY')

    def fill_trigger_id_unphysical(self, value=-1):
        """Fill the TRG_ID column with an unphysical value.

        This is called once when the event lists for the single source components
        are simulated in order to have the column in the event list itself
        filles and of the proper length. The actual trigger IDs are then set to
        the actual values after the source components are merged into the final
        event list and the deadtime is applied.
        """
        self._set_column('TRG_ID', value)

    def fill_trigger_id(self):
        """Fill the TRG_ID column with sequential values starting from 1.
        """
        self._set_column('TRG_ID', numpy.arange(1, self.num_events() + 1))

    def fill_livetime_unphysical(self, value=-1):
        """Fill the LIVETIME column with an unphysical value.

        This is called once when the event lists for the single source components
        are simulated in order to have the column in the event list itself
        filles and of the proper length. The actual trigger IDs are then set to
        the actual values after the source components are merged into the final
        event list and the deadtime is applied.
        """
        self._set_column('LIVETIME', value)

    def fill_livetime(self, gti_list, deadtime=0.):
        """Fill the LIVETIME column.

        The basic rule is that the livetime for the first event is the time
        elapsed since the start run, and for all the other events we just calculate
        the time difference and subtract the deadtime per event.

        Args
        ----
        gti_list : xGTIList instance
            The list of good time intervals for the observation.

        deadtime : float
            The deadtime per event.
        """
        # Recover the start time of the observation and of all the good
        # time intervals (start_met will be equal to the first element of
        # start_mets).
        start_met = gti_list.start_met
        start_mets = numpy.array(gti_list.start_mets())
        # Cache the list of event times. Note that we prepend the start time of
        # the observation so that the following numpy.diff() has the proper length.
        met = numpy.append(start_met, self.time())
        # Find the indices of the first events in each good time interval.
        # Note we bisect on the right, so that the first index will be 1.
        idx = numpy.searchsorted(met, start_mets, side='right')
        # Calculate the time differences between successive events.
        dt = numpy.diff(met)
        # Overwrite the time difference for the first event in each time interval.
        # The extra +deadtime at the end is such that when we the subtract the
        # deadtime at the following line, this cancels out for this small subset
        # of events.
        dt[idx - 1] = met[idx] - start_mets + deadtime
        # Subtract the deadtime.
        dt -= deadtime
        # Round up and convert to us.
        dt = numpy.floor(dt * 1.e6)
        self._set_column('LIVETIME', dt)

    def set_num_clusters(self, value=1):
        """Fill the NUM_CLU column with a fixed value.

        Arguments
        ---------
        value : number or array
            The values for the column to be filled with.

        Note the control on whether value is a number or an array is performed in
        the _set_column() method.
        """
        self._set_column('NUM_CLU', value)

    def set_weights(self, value=1.):
        """Fill the W_MOM column.

        Arguments
        ---------
        value : number or array
            The values for the column to be filled with.

        Note the control on whether value is a number or an array is performed in
        the _set_column() method.
        """
        self._set_column('W_MOM', value)

    def set_mc_energy_columns(self, mc_energy, mc_pha, mc_pi):
        """Set all the columns related to the Monte Carlo energy.
        """
        self._set_column('MC_ENERGY', mc_energy)
        self._set_column('MC_PHA', mc_pha)
        self._set_column('MC_PI', mc_pi)

    def set_energy_columns(self, energy, pha, pi):
        """Set all the columns related to the measured energy.
        """
        self._set_column('ENERGY', energy)
        self._set_column('PHA', pha)
        self._set_column('PI', pi)

    def set_mc_sky_position_columns(self, mc_ra, mc_dec, mc_x, mc_y):
        """Set all the columns related to the Monte Carlo sky position.
        """
        self._set_column('MC_RA', mc_ra)
        self._set_column('MC_DEC', mc_dec)
        self._set_column('MC_X', mc_x)
        self._set_column('MC_Y', mc_y)

    def set_sky_position_columns(self, ra, dec, x, y):
        """Set all the columns related to the measured sky position.
        """
        self._set_column('RA', ra)
        self._set_column('DEC', dec)
        self._set_column('X', x)
        self._set_column('Y', y)

    def set_detector_position_columns(self, detx, dety):
        """Set all the columns related to the event position in the GPD frame.
        """
        self._set_column('DETX', detx)
        self._set_column('DETY', dety)

    def set_phi_columns(self, phi, detphi):
        """Set all the columns related to the photoelectron direction.
        """
        self._set_column('PHI', phi)
        self._set_column('DETPHI', detphi)
        self._set_column('Q', xStokesAnalysis.stokes_q(phi, weights=None))
        self._set_column('U', xStokesAnalysis.stokes_u(phi, weights=None))

    def set_mc_gain_column(self, gain=1.):
        """Set the column corresponding to the relative GEM gain.
        """
        self._set_column('MC_GAIN', gain)

    def set_seed_columns(self, mc_energy, mc_pha, mc_pi, mc_ra, mc_dec,
                         mc_x, mc_y, phi, detphi):
        """Set all the columns for the seed event list.

        These is the first set of columns being simulated in xpobssim, which is
        why we have this convenience function to set them all at once.

        Note
        ----
        Historically the columns set by this function included the time, but
        that was abandoned for the sake of consistency when we added the
        possibility of setting the event times at construction time. So the
        basic workflow is now:

        * initialize the event list using right after the event times have been
          extracted passing the times themselves and the source identifier to
          the constructor;
        * call set_seed_columns();
        * call set_rec_columns().
        """
        self.set_mc_energy_columns(mc_energy, mc_pha, mc_pi)
        self.set_mc_sky_position_columns(mc_ra, mc_dec, mc_x, mc_y)
        self.set_phi_columns(phi, detphi)

    def set_rec_columns(self, pha, pi, energy, ra, dec, x, y, detx, dety):
        """Set all the columns of the EVENT extension that imply a
        convolution of the MC truth with the instrument response functions.
        """
        self.set_energy_columns(energy, pha, pi)
        self.set_sky_position_columns(ra, dec, x, y)
        self.set_detector_position_columns(detx, dety)

    def apply_fiducial_area(self):
        """Trim out the events falling outside the detector active area.
        """
        detx, dety = self.detector_coordinates()
        self.apply_fiducial_area_base(detx, dety)

    def apply_vignetting(self, vign, ra_pnt, dec_pnt):
        """Apply the effective area vignetting to the event list.
        """
        ra, dec = self.mc_sky_coordinates()
        energy = self.mc_energy()
        self.apply_vignetting_base(ra, dec, energy, vign, ra_pnt, dec_pnt)

    @timing
    def apply_charging(self, irf_set, **kwargs):
        """Apply the GEM charging to the event list.
        """
        logger.info('Applying GEM charging to the event list...')
        # Read the initial charging map file
        try:
            charging_map_file = kwargs.get('chrgmaps')[irf_set.du_id]
        except KeyError:
            abort('No matching charging map file for du %s' % irf_set.du_id)
        initial_dg_fast, initial_dg_slow = read_charging_map(charging_map_file)
        # Create the energy-rate cube object to perform the correction.
        nside = initial_dg_fast.shape[0]
        start_met = kwargs.get('start_met')
        stop_met = kwargs.get('stop_met')
        met_step = kwargs.get('chrgtstep')
        assert met_step > 0
        num_met_steps = round((stop_met - start_met)/met_step)
        if num_met_steps <= 1:
            met_edges = numpy.array([start_met, stop_met])
        else:
            # We define our time edges so that they are equally spaced and
            # end exactly at 'tstop', even if that means that their
            # size is not exactly equal to 'met_step'
            met_edges = numpy.linspace(start_met, stop_met, int(num_met_steps))
        cube = xEnergyFluxCube(nside, met_edges)
        # Fill the cube---mind that, since we are doing this before throwing
        # away the events due to the dead time, we are setting the
        # deadtime correction to None. The other thing that might be worth
        # noticing is whether we should fill the cube with the true or measured
        # energies (not that it will make a huge difference). Since the
        # measured energy is representative of the charge fluctuation in the
        # primary ionization and subsequent avalache, I'd argue that this is
        # *exactly* what determines the charging and should be used here.
        detx, dety = self.detector_coordinates()
        time_ = self.time()
        energy = self.energy()
        cube.fill(detx, dety, time_, energy, deadtime_correction=None)
        # And, finally, modify the event energies and related quantities.
        # Read the charging parameters file
        try:
            charging_params_file = kwargs.get('chrgparams')[irf_set.du_id]
        except KeyError:
            abort('No matching charging parameter file for du %s' \
                 % irf_set.du_id)
        k_c_fast, tau_d_fast, delta_max_fast, k_c_slow, tau_d_slow, delta_max_slow = \
            read_charging_parameters(charging_params_file)
        cube.calculate_gain_data(k_c_fast, tau_d_fast, delta_max_fast, initial_dg_fast,
            k_c_slow, tau_d_slow, delta_max_slow, initial_dg_slow)
        gain = cube.gain(detx, dety, time_)
        self.set_mc_gain_column(gain)
        # Now we are ready to change all the variables related to the measured
        # energy (energy, pha and pi).
        corr_energy = energy * gain
        corr_pha, corr_pi = irf_set.edisp.pha_analysis(corr_energy)
        self.set_energy_columns(corr_energy, corr_pha, corr_pi)
        logger.info('Done, average energy scaling: %.3f', gain.mean())
        self.__charging_maps = cube.charging_map

    @timing
    def apply_dead_time(self, deadtime):
        """Apply the dead time effect to the event list.
        """
        logger.info('Applying dead time (%.5f) to the event list...', deadtime)
        num_events = self.num_events()
        evt_time = self.time()
        good = numpy.ones(evt_time.shape, dtype=bool)
        last_time = evt_time[0]
        for i in range(1, num_events):
            if (evt_time[i] - last_time) < deadtime:
                good[i] = False
            else:
                last_time = evt_time[i]
        self.trim(good)
        frac = 100. * self.num_events() / float(num_events)
        logger.info('Done, %d/%d (%.2f%%) events remaining', self.num_events(), num_events, frac)

    def _finalize(self, irf_set, **kwargs):
        """Finalize the event list.

        This should be run when there are no more events to be added to the list
        and we're ready to write the list itself to the output file. The method
        basically runs several different task, if and when necessary---in this order:

        * run `apply_fiducial_area()`
        * sort the event list;
        * run `apply_charging()`;
        * run `apply_dead_time()`;
        * run `fill_livetime()`
        * run `fill_trigger_id()`.

        Note that, at simulation time, we do want to do calculate the charging
        *after* the event list is complete, the vignetting and the fiducial
        area cuts (if any) have been applied, but *before* the deadtime is
        applied, so that we don't have to bother about any deadtime correction
        (the charging cares about the real rate of events impinging onto the
        detector, not only those that actually trigger).

        On the other hand, the LIVETIME and TRIGGER_ID columns should be filled
        at the every end, when the final list of events that have been read out
        is frozen.
        """
        # If there are no events in the list we refrain from doing anything.
        if self.num_events() == 0:
            return
        logger.info('Finalizing the event list...')
        self.apply_fiducial_area()
        self.sort()
        if kwargs.get('charging'):
            self.apply_charging(irf_set, **kwargs)
        deadtime = kwargs.get('deadtime')
        if deadtime > 0.:
            self.apply_dead_time(deadtime)
        gti_list = kwargs.get('gti_list')
        self.fill_livetime(gti_list, deadtime)
        self.fill_trigger_id()

    def write_fits(self, creator, roi_model, irf_set, **kwargs):
        """Write the event list and associated ancillary information to file.
        """
        du_id = irf_set.du_id
        self._finalize(irf_set, **kwargs)
        file_path = kwargs.get('outfile')
        irf_name = kwargs.get('irfname')
        start_met = kwargs.get('start_met')
        duration = kwargs.get('duration')
        stop_met = start_met + duration
        gti_list = kwargs['gti_list']
        ontime = gti_list.total_good_time()
        livetime = 1.e-6 * self.livetime()
        # If the on-orbit calibration is enabled we need to limit the sum of the
        # livetimes to the events from the celestial source.
        if kwargs.get('onorbitcalib'):
            ocsrcid = kwargs.get('calib_pattern').calibration_source.identifier
            livetime = livetime[self.src_id() != ocsrcid]
        livetime_sum = livetime.sum()
        deadtime_correction = livetime_sum / ontime

        # Cache the proper args for the common header keywords.
        _time_args = start_met, stop_met, duration, ontime, livetime_sum, deadtime_correction
        _obj_args = roi_model.ra, roi_model.dec, kwargs.get('objname')

        def _update_header(hdu):
            """Small nested functions to update all the header keywords that are
            relevant for all the extensions.
            """
            set_telescope_header_keywords(hdu, du_id)
            set_time_header_keywords(hdu, *_time_args)
            set_object_header_keywords(hdu, *_obj_args)

        # Create the primary header.
        primary_hdu = xLvl2PrimaryHDU(creator=creator)
        _update_header(primary_hdu)

        # Create the EVENTS extension.
        event_hdu = xBinTableHDUEvents(roi_model.ra, roi_model.dec, self)
        _update_header(event_hdu)
        # And, as we discovered in issue https://github.com/lucabaldini/ixpeobssim/issues/688
        # we do need to set TLMIN and TLMAX for the PI column.
        set_tlbounds(event_hdu, 'PI', TLMIN, TLMAX)
        # Support for pseudo-Lv1a files for interoperability with the official
        # IXPE charging correction tool, see
        # https://bitbucket.org/ixpesw/ixpeobssim/issues/392
        if kwargs.get('lv1a'):
            # Tweak the header keywords.
            version = kwargs.get('lv1version')
            set_version_keywords(primary_hdu, version)
            set_version_keywords(event_hdu, version)
            # Calculate the PHA_EQ column. We arbitrarily set the PHA_EQ to
            # 18000 at the energy of 5.9 keV.
            pha_eq = self.energy() * 18000. / 5.9
            pha_eq = fits.Column(name='PHA_EQ', format='D', array=pha_eq)
            detx, dety = self.detector_coordinates()
            barx = fits.Column(name='BARX', format='D', array=detx)
            bary = fits.Column(name='BARY', format='D', array=dety)
            additional_cols = fits.ColDefs([barx, bary, pha_eq])
            cols = event_hdu.columns + additional_cols
            event_hdu = fits.BinTableHDU.from_columns(cols, header=event_hdu.header)

        # Create the MONTE_CARLO extension.
        mc_hdu = xBinTableHDUMonteCarlo(self)
        mc_hdu.set_irf_name(irf_name)
        _update_header(mc_hdu)

        # Create the GTI extension.
        gti_hdu = xBinTableHDUGTI([gti_list.start_mets(), gti_list.stop_mets()])
        _update_header(gti_hdu)

        # Create the ROITABLE extension
        _src_id = numpy.array([src.identifier for src in roi_model.values()])
        _src_name = numpy.array([src.name for src in roi_model.values()])
        roi_hdu = xBinTableHDURoiTable([_src_id, _src_name])
        roi_hdu.set_center(roi_model.ra, roi_model.dec)
        _update_header(roi_hdu)

        # Preparing the HDU list with the mandatory extensions.
        hdu_list = fits.HDUList([primary_hdu, event_hdu, mc_hdu, gti_hdu, roi_hdu])

        # Create the TIMELINE extension.
        if kwargs.get('timelinedata'):
            logger.info('Creating the TIMELINE extension...')
            epoch_data = kwargs.get('timeline').epoch_data()
            timeline_hdu = xBinTableHDUTimeline(epoch_data)
            _update_header(timeline_hdu)
            hdu_list.append(timeline_hdu)

        # Create the SC_DATA extension.
        if kwargs.get('scdata'):
            logger.info('Creating the SC_DATA extension...')
            dither_params = parse_dithering_kwargs(**kwargs)
            args = kwargs.get('scdatainterval'), dither_params, kwargs.get('saa'), \
                kwargs.get('occult')
            sc_data = kwargs.get('timeline').sc_data(*args)
            scdata_hdu = xBinTableHDUSpacecraftData(sc_data)
            scdata_hdu.set_roll_angle(kwargs.get('roll'))
            _update_header(scdata_hdu)
            hdu_list.append(scdata_hdu)

        # Create the OCTI extension.
        if kwargs.get('onorbitcalib'):
            logger.info('Creating the OCTI extension...')
            calib_pattern = kwargs.get('calib_pattern')
            octi_data = calib_pattern.octi_data(du_id)
            octi_hdu = xBinTableHDUOCTI(octi_data)
            octi_hdu.set_cal_stats(calib_pattern.num_calibration_runs(du_id),
                calib_pattern.total_calibration_time(du_id))
            _update_header(octi_hdu)
            hdu_list.append(octi_hdu)

        # Add a CHRG_MAP extesion
        if kwargs.get('charging'):
            logger.info('Creating the CHRG_MAP extension...')
            map_date = met_to_string_utc(stop_met, fmt="%m/%d/%Y")
            map_time = met_to_string_utc(stop_met, fmt="%H:%M:%S")
            fast_map = self.__charging_maps[0]
            slow_map = self.__charging_maps[1]
            _kwargs = dict(start_date=map_date, start_time=map_time)
            charging_hdu = create_charging_map_extension(fast_map, slow_map, **_kwargs)
            hdu_list.append(charging_hdu)

        # We're good to go.
        hdu_list.info()
        hdu_list.writeto(file_path, overwrite=True)
        hdu_list.close()
        logger.info('Event list written to %s...', file_path)



class xEventFile:

    """Read/write interface to level-2 event files.

    This is the main class providing access to the information into a
    FITS file containing an event list.

    It supports minimal output facilities, such as adding columns and
    writing to file.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        check_input_file(file_path, 'fits')
        logger.info('Opening input event file %s...', file_path)
        self.hdu_list = fits.open(file_path)
        self.hdu_list.info()
        self.primary_header = self.hdu_list['PRIMARY'].header
        self.event_data = self.hdu_list['EVENTS'].data
        self._wcs = self._retrieve_wcs_info()
        logger.info(self._wcs)
        self.mc_data = self._retrieve_mc_data()
        self.roi_table = self._build_roi_table()

    def _retrieve_wcs_info(self):
        """Read the WCS information from the EVENTS header.

        .. warning::

           Need to fully understand why we have to manually set the array_shape
           property of the WCS read from the EVENTS header.
        """
        logger.info('Reading WCS information from the EVENTS HDU...')
        try:
            wcs_ = wcs.WCS(self.hdu_list['EVENTS'].header, keysel=['pixel'])
            wcs_.array_shape = (_SKYCOORD_NUM_SIDE_PIXELS, _SKYCOORD_NUM_SIDE_PIXELS)
            return wcs_
        except Exception as e:
            logger.error(e)
            logger.error('Cannot parse the WCS information.')
            logger.warning('Creating a standard WCS on the fly based on RA_OBJ and DEC_OBJ...')
            logger.warning('This might get you going, but you should fix the input file.')
            ra0 = self.primary_header.get('RA_OBJ')
            dec0 = self.primary_header.get('DEC_OBJ')
            return build_standard_wcs(ra0, dec0)

    def _retrieve_mc_data(self):
        """Read the MONTE_CARLO extension, if available.
        """
        logger.info('Loading data from the MONTE_CARLO extension...')
        try:
            return self.hdu_list['MONTE_CARLO'].data
        except KeyError:
            logger.info('Event file comes without a MONTE_CARLO extension.')
            return None

    def _build_roi_table(self):
        """Rebuild the ROI table based in the information in the ROITABLE
        extension of the event file.
        """
        logger.info('Re-building the ROI table...')
        try:
            data = self.hdu_list['ROITABLE'].data
            return {id_: name for id_, name in zip(data['SRCID'], data['SRCNAME'])}
        except KeyError:
            logger.info('Input file has no ROITABLE extension, cannot rebuild ROI table.')
            return {}

    def close(self):
        """Close the underlying HDU.
        """
        logger.info('Closing event file %s...', self.file_path())
        self.hdu_list.close()

    def __del__(self):
        """Close the underlying HDU upon object destruction.

        It is not entirely clear to me whether this is really needed, but I do
        see some ResourceWarning in the unit tests without this custom
        destructor. This might warrant some detailed study.
        """
        self.close()

    def file_path(self):
        """Return the path to the underlying file.
        """
        return self.hdu_list.filename()

    def num_events(self):
        """Return the total number of events in the event file.
        """
        return len(self.event_data)

    def wcs_reference(self):
        """Return the reference point of the underlying WCS object.

        Note that, since we are making sure that an event file *always* has a
        valid WCS associated with it (be it the actual content of the EVENTS
        header or one created on the fly given RA_SRC and DEC_SRC), this
        should be the preferred way to gauge the center of the field---compared,
        e.g., to using RA_SRC, DEC_SRC or RA_PNT, DEC_PNT directly.
        """
        return tuple(self._wcs.wcs.crval)

    def primary_keywords(self):
        """Return a list of all the relevant keywords in the PRIMARY header.

        This is used, e.g., to propagate all the necessary information (such as
        the ROI and the IRFs used) from the event files to the binned
        data files that are created from them.

        Warning
        -------
        This function needs to be refactored!
        """
        _header = self.primary_header
        keywords = []
        for key in ['TELESCOP', 'INSTRUME', 'DETNAM'] + \
            [item[0] for item in _TIME_HEADER_KEYWORDS]:
            keywords.append((key, _header[key], _header.comments[key]))
        keywords.append(('DETCHANS', NUM_CHANNELS, ''))
        backscal = self.backscal()
        if backscal is not None:
            keywords.append(('BACKSCAL', backscal, ''))
        return keywords

    def start_met(self):
        """Return the start MET of the observation.
        """
        return self.primary_header['TSTART']

    def stop_met(self):
        """Return the stop MET of the observation.
        """
        return self.primary_header['TSTOP']

    def min_good_time(self):
        """Return the smaller START time for the GTIs in the event file.
        """
        return self.hdu_list['GTI'].data['START'].min()

    def max_good_time(self):
        """Return the largest STOP time for the GTIs in the event file.
        """
        return self.hdu_list['GTI'].data['STOP'].max()

    def total_good_time(self):
        """Return the sum of all the GTIs in the event file.
        """
        start = self.hdu_list['GTI'].data['START']
        stop = self.hdu_list['GTI'].data['STOP']
        return (stop - start).sum()

    def livetime(self):
        """Return the livetime, i.e. the duration corrected for dead time.
        """
        return self.primary_header['LIVETIME']

    def deadtime_correction(self):
        """Return the deadtime correction.
        """
        return self.primary_header['DEADC']

    def backscal(self):
        """Return the value of the BACKSCAL header keyword, if present.

        The keyword is written out by xpselect when appropriate. The function
        returns None when the keyword is missing in the event file.
        """
        try:
            return self.primary_header['BACKSCAL']
        except KeyError:
            return None

    def irf_name(self):
        """Return the name of the IRF set used to run the simulation.

        Note the try/except block in the function body was added for ixpeobssim
        to inter-operate with the ixpesim output photon lists, which might
        have a MONTE_CARLO extension without the IRFNAME header keyword.
        In this case you have to set the keyword by hand, e.g., passing the
        proper command-line switch to xpbin.
        """
        if self.mc_data is None:
            msg = 'MONTE_CARLO extension not available, cannot retrieve simulation IRF name'
            logger.warning(msg)
            return None
        try:
            return self.hdu_list['MONTE_CARLO'].header['IRFNAME']
        except KeyError:
            logger.warning('MONTE_CARLO extension has no "IRFNAME" keyword')
            return None

    def du_id(self):
        """Return the detector unit number used to run the simulation.
        """
        return int(self.primary_header['DETNAM'][-1])

    def time_data(self):
        """Return the TIME column.
        """
        return self.event_data['TIME']

    def livetime_data(self):
        """Return the LIVETIME column.
        """
        return self.event_data['LIVETIME']

    def phase_data(self):
        """Return the PHASE column.
        """
        try:
            return self.event_data['PHASE']
        except KeyError:
            logger.warning('Event file has no PHASE column, you might need to run xpphase.py...')
            return None

    def pi_data(self, mc=False):
        """Return the PI column.
        """
        if not mc:
            return self.event_data['PI']
        if self.mc_data is None:
            logger.warning('MONTE_CARLO extension not available, cannot retrieve MC_PI data...')
            return None
        return self.mc_data['MC_PI']

    def energy_data(self, mc=False):
        """Return the ra and dec column data, in either the Monte Carlo or
        the measured flavor depending on the keyword arguments passed to the
        binning class.

        The try/except clause if needed for inter-operability with the ixpesim
        output files, where the column naming in the MONTE_CARLO extension is
        not in line with ixpeobssim, see
        https://bitbucket.org/ixpesw/gpdsw/issues/327
        """
        if not mc:
            return channel_to_energy(self.pi_data())
        if self.mc_data is None:
            logger.warning('MONTE_CARLO extension not available, cannot retrieve MC_ENERGY data...')
            return None
        try:
            return self.mc_data['MC_ENERGY']
        except KeyError:
            return self.mc_data['ENERGY']

    def xy_data(self):
        """Return the X and Y columns.
        """
        return self.event_data['X'], self.event_data['Y']

    def sky_position_data(self, mc=False):
        """Return the ra and dec column data, in either the Monte Carlo or
        the measured flavor depending on the keyword arguments passed to the
        binning class.

        Note that, since flight level-2 data do not contain RA and DEC columns,
        this uses internally X and Y.
        """
        if not mc:
            return self._wcs.wcs_pix2world(*self.xy_data(), WCS_ORIGIN)
        if self.mc_data is None:
            abort('MONTE_CARLO extension not available, cannot retrieve MC_RA/DEC data...')
        try:
            return self.mc_data['MC_RA'], self.mc_data['MC_DEC']
        except KeyError:
            abort('MC_RA and/or MC_DEC columns not available...')

    def q_data(self):
        """Return the Q column.
        """
        return self.event_data['Q']

    def u_data(self):
        """Return the U column.
        """
        return self.event_data['U']

    def stokes_data(self):
        """Return the Q and U columns.
        """
        return self.q_data(), self.u_data()

    def det_position_data(self):
        """Return the detx and dety column data, in either the detector x and y
        measured coordinates on the detector.
        """
        return self.event_data['DETX'], self.event_data['DETY']

    def phi_data(self):
        """Return the PHI column.
        """
        return self.event_data['PHI']

    def trigger_id_data(self):
        """Return the TRG_ID column.
        """
        return self.event_data['TRG_ID']

    def srcid_data(self):
        """Return the SRC_ID column.

        Note the try/except block in the function body was added for ixpeobssim
        to inter-operate with the ixpesim output photon lists, which come
        have a MONTE_CARLO extension without the ``SCR_ID`` column.
        """
        if self.mc_data is None:
            logger.warning('MONTE_CARLO extension not available, cannot retrieve SRC_ID data')
            return None
        try:
            return self.mc_data['SRC_ID']
        except KeyError:
            logger.warning('MONTE_CARLO extension has no SRC_ID column')
            return None

    def ds9_region_mask(self, *region_list, mc=False):
        """Check which events are inside the logical or of all the regions
        within a given region list and return the corresponding array mask.
        """
        ra, dec = self.sky_position_data(mc)
        return ds9_region_filter_sky(ra, dec, self._wcs, *region_list)

    def ds9_region_file_mask(self, file_path, mc=False):
        """Convenience function filtering a photon list from a ds9 region file.
        """
        logger.info('Loading regions from %s...', file_path)
        region_list = read_ds9(file_path)
        return self.ds9_region_mask(*region_list, mc=mc)

    def _extension_data(self, ext_name):
        """Convenience function to retrieve the extension data.
        """
        try:
            data = self.hdu_list[ext_name].data
        except KeyError:
            logger.warning('Input file has no extension %s', ext_name)
            return None
        return {col_name: data[col_name] for col_name in data.columns.names}

    def timeline_data(self):
        """Return the timeline data.
        """
        return self._extension_data(xBinTableHDUTimeline.NAME)

    def sc_data(self):
        """Return the spacecraft data.
        """
        return self._extension_data(xBinTableHDUSpacecraftData.NAME)

    def gti_data(self):
        """Return thr GTI data.
        """
        return self._extension_data(xBinTableHDUGTI.NAME)

    def octi_data(self):
        """Return the OCTI data.
        """
        return self._extension_data(xBinTableHDUOCTI.NAME)

    def average_deadtime_per_event(self):
        """Calculate the average deadtime per event.

        .. warning::

           Note this is trickier than it seems, as the result might not be accurate
           if the underlying event file has already been filteres in time and/or phase.
        """
        # Calculate the average deadtime per event for the full file. Note
        dead_time = self.primary_header.get('ONTIME') - self.primary_header.get('LIVETIME')
        logger.info('Total dead time: %.3f s', dead_time)
        dtpe = dead_time / self.num_events()
        logger.info('Estimated average dead time per event: %.3f ms', 1000. * dtpe)
        # Minimal cross-ckeck on the deadtime calculation---we want it to be
        # *at the very least* larger than the minimum delta event time.
        min_delta_time = numpy.diff(self.time_data()).min()
        logger.info('Minimum delta event time: %.3f ms', 1000. * min_delta_time)
        if dtpe < min_delta_time:
            logger.warning('Average deadtime smaller than the smallest delta event time...')
            logger.warning('This might indicate that the original file has been filtered!')
            logger.warning('You might want to double check the livetime calculation manually :-)')
        return dtpe

    def pol_deg_weighted_average(self, pol_deg_model, ebins=10):
        """Calculate the weighted average of the polarization degree as a
        function of energy on the given event sample for a generic input
        polarization model.

        This is achieved by calculating the input model of the polarization
        degree on an event by event basis (i.e., evaluating the model itself on
        the proper columns of the event file), creating a 2-dimensional
        histogram of the polarization degree vs. energy, and calculating the
        average of the polarization degree in vertical slices of energy.

        Warning
        -------
        The polarization degree is averaged coherently, i.e., under the
        assumption that the polarization angle is constant over all the event
        sample. Therefore using this function only makes sense when either the
        input polarization model has a constant polarization angle, or after the
        photoelectron directions in the event file have been properly aligned to
        the given input model.

        Args
        ----
        pol_deg_model : callable
            A Python function or an otherwise callable object with the proper
            signature returning the polarization degree as a function of
            energy, time and sky position (in this order).

        ebins : int or arra-like (default None)
            The energy binning for the weighted average. Loosely modeled on the
            numpy and scipy digitization functions, if bins is and integer,
            a logarithmically spaced binning between 2 and 8 keV is created.

        Return
        ------
        An xScatterPlot object encoding the average polarization degree as a
        function of energy.
        """
        # Retrieve the necessary data columns (energy, time, and sky position)
        energy = self.energy_data()
        time_ = self.time_data()
        ra, dec = self.sky_position_data()
        # Calculate the input model on the event sample.
        pol_deg = pol_deg_model(energy, time_, ra, dec)
        # Create the temporary 2-dimensional histogram.
        if isinstance(ebins, int):
            ebins = numpy.logspace(numpy.log10(2.), numpy.log10(8.), ebins)
        assert isinstance(ebins, numpy.ndarray)
        hist_2d = xHistogram2d(ebins, numpy.linspace(0., 1., 250))
        hist_2d.fill(energy, pol_deg)
        # Project the vertical columns and return the scatter plot with the
        # average polarization degree values.
        x = hist_2d.bin_centers(0)
        y = numpy.array([h.mean(0) for h in hist_2d.vslices()])
        return xScatterPlot(x, y, xlabel='Energy [keV]', ylabel='Polarization degree')

    def filter(self, mask, history=None):
        """Filter the event file according to a given mask of selected rows.

        This filtering happens in place, so it should be used when you don't
        need to apply subsequent selections on the same file multiple times.
        """
        logger.info('Filtering in place event file (%d out of %d rows selected)',
                    mask.sum(), len(mask))
        primary_header = self.hdu_list['PRIMARY'].header
        primary_header['DATE'] = current_datetime_string_utc()
        if history is not None:
            primary_header['HISTORY'] = history
        self.hdu_list['EVENTS'].data = self.hdu_list['EVENTS'].data[mask]
        if self.mc_data is not None:
            self.hdu_list['MONTE_CARLO'].data = self.hdu_list['MONTE_CARLO'].data[mask]
        return self.hdu_list

    def copy_and_filter(self, mask, history=None):
        """Filter the event file according to a given mask of selected rows.

        This is similar in spirit to the filter() method, except that the
        filtering is not done in place, but everything is copied beforehand,
        so that the original HDU list is not modified, and multiple selections
        can be applied in series, at the expense of a larger memory footprint.
        """
        logger.info('Copy-filtering event file (%d out of %d rows selected)',
                    mask.sum(), len(mask))
        primary_header = self.hdu_list['PRIMARY'].header.copy()
        primary_header['DATE'] = current_datetime_string_utc()
        if history is not None:
            primary_header['HISTORY'] = history
        hdu_list = [fits.PrimaryHDU(header=primary_header)]
        hdu = self.hdu_list['EVENTS'].copy()
        hdu.data = hdu.data[mask]
        hdu_list.append(hdu)
        if self.mc_data is not None:
            hdu = self.hdu_list['MONTE_CARLO'].copy()
            hdu.data = hdu.data[mask]
            hdu_list.append(hdu)
        if self.roi_table:
            hdu_list.append(self.hdu_list['ROITABLE'].copy())
        return fits.HDUList(hdu_list)

    def set_column(self, ext_name, col_name, col_data):
        """Overwrite the data for an existing column.

        Args
        ----
        ext_name : str
            The name of the extension that the column should be added to.

        col_name : str
            The name of the column to be added.

        col_format : str
            The format for the new column (e.g., 'E').

        col_data : array_like
            The data for the new column.
        """
        logger.info('Overwriting column %s for the event file...', col_name)
        self.hdu_list[ext_name].data[col_name] = col_data

    def add_column(self, ext_name, col_name, col_format, col_data):
        """Add a column to the specified extension of the event file.

        Args
        ----
        ext_name : str
            The name of the extension that the column should be added to.

        col_name : str
            The name of the column to be added.

        col_format : str
            The format for the new column (e.g., 'E').

        col_data : array_like
            The data for the new column.
        """
        logger.info('Adding column %s to the event file...', col_name)
        hdu = self.hdu_list[ext_name]
        cols = hdu.data.columns
        cols += fits.Column(name=col_name, format=col_format, array=col_data)
        hdu = fits.BinTableHDU.from_columns(cols, header=hdu.header)
        self.hdu_list[ext_name] = hdu

    def add_columns(self, ext_name, *columns):
        """Add multiple columns to the specified extension of the event file.

        Compared to multiple calls to add_column(), and at the cost of a slight
        code duplication, this avoids create a new binary table at each call.

        Args
        ----
        ext_name : str
            The name of the extension that the column should be added to.

        columns : fits.Column instances
            The columns to be added.
        """
        logger.info('Adding columns %s to the event file...', [col.name for col in columns])
        hdu = self.hdu_list[ext_name]
        cols = hdu.data.columns + columns
        hdu = fits.BinTableHDU.from_columns(cols, header=hdu.header)
        self.hdu_list[ext_name] = hdu

    def remove_columns(self, ext_name, *col_names):
        """Remove one or more columns from the specified extension of the event file.

        Args
        ----
        ext_name : str
            The name of the extension that the column should be added to.

        *col_names : str
            The name(s)  of the column(s) to be removed
        """
        logger.info('Removing columns %s from the event file...', col_names)
        hdu = self.hdu_list[ext_name]
        cols = [col for col in hdu.data.columns if col.name not in col_names]
        hdu = fits.BinTableHDU.from_columns(cols, header=hdu.header)
        self.hdu_list[ext_name] = hdu

    def write(self, file_path, overwrite=False):
        """Write the underlying HDU list to file.
        """
        logger.info('Writing event data %s...', file_path)
        self.hdu_list.info()
        self.hdu_list.writeto(file_path, overwrite=overwrite)
        logger.info('Done.')

    def write_fits_selected(self, selection_mask, file_path, header_keywords=None,
                            history=None, overwrite=True, filter_in_place=True):
        """Write to file a subselection of events.

        Arguments
        ---------

        selection_mask : array
            A mask for the row selection.

        file_path : string
            The path to the output file.
        """
        if filter_in_place:
            hdu_list = self.filter(selection_mask, history)
        else:
            hdu_list = self.copy_and_filter(selection_mask, history)
        if header_keywords is not None:
            for hdu in hdu_list:
                if hdu.name in ('PRIMARY', 'EVENTS', 'GTI'):
                    logger.debug('Updating keywords in %s extension', hdu.name)
                    for key, value in header_keywords.items():
                        hdu.header.set(key, value)
        logger.info('Writing data subselection to %s...', file_path)
        hdu_list.info()
        hdu_list.writeto(file_path, overwrite=overwrite)
        logger.info('Done.')


class xEventFileFriend:

    """Read/write interface to put together level-1 and level-2 event files.

    It supports minimal output facilities, mainly re-implementing functions
    from xEventFile class
    """

    def __init__(self, file_list2, file_list1=None):
        """Simple wrapper class to read two list of files, Lv2 and Lv1,
        for one observation ad get variables for only events filtered
        in the level2 selection
        """
        if not isinstance(file_list2, list):
            file_list2 = [file_list2]
        if not isinstance(file_list1, list):
            file_list1 = [file_list1]
        self.file_list2 = []
        self.file_list1 = []
        self.time_ids = None
        time1 = numpy.array([])
        time2 = numpy.array([])
        for f2 in file_list2:
            self.file_list2.append(xEventFile(f2))
            time2 = numpy.append(time2, self.file_list2[-1].time_data())
        if None in file_list1:
            # just concatenate Level2  - mainly for MC use
            self.file_list1 = None
        else:
            for f1 in file_list1:
                self.file_list1.append(xEventFile(f1))
                time1 = numpy.append(time1, self.file_list1[-1].time_data())
            # create a mask to select Lv1 on for the whole lists of files
            # True if Lv1 evt time is also in Lv2 evt time list
            # assume Lv1 list is larger than Lv2
            #self.time_ids.append(numpy.nonzero(numpy.in1d(time1, time2))[0])
            self.time_ids = numpy.in1d(time1, time2)

    def l1value(self, val, all_events=False):
        """
        """
        if self.file_list1 == None:
            return None
        outvalues = numpy.array([])
        for fl1 in self.file_list1:
            outvalues = numpy.append(outvalues, fl1.event_data[val])
        if all_events:
            return outvalues
        else:
            return outvalues[self.time_ids]

    def l2value(self, val):
        """
        """
        values = None
        for fl2 in self.file_list2:
            if values is None:
                values =  fl2.event_data[val]
            else:
                values =  numpy.append(values, fl2.event_data[val])
        return values

    def energy_data(self, mc=False):
        """Wrap xEventFile.energy_data() appending values for all the LV2 files
        """
        values = numpy.array([])
        for fl2 in self.file_list2:
            values = numpy.append(values, fl2.energy_data(mc))
        return values

    def energy_l1_data(self):
        """Rewrite xEventFile.energy_data() to get energy in keV for L1 data
        """
        logger.warning('xEventFileFriend.energy_l1_data() not implemented yet')
        return None
        # if self.file_list1 == None:
        #     return None
        # outvalues = None
        # for (i, fl1) in enumerate(self.file_list1):
        #     values = channel_to_energy(fl1.event_data['PI']) # NO! need conversion from PHA
        #     if outvalues is None:
        #         outvalues =  values
        #     else:
        #         outvalues =  numpy.append(outvalues, values)
        # return outvalues

    def sky_position_data(self, mc=False):
        """Wrap xEventFile.sky_position_data() appending values for all the LV2 files
        """
        ra = numpy.array([])
        dec = numpy.array([])
        for fl2 in self.file_list2:
            val1, val2 = fl2.sky_position_data(mc)
            ra  = numpy.append(ra, val1)
            dec = numpy.append(dec, val2)
        return ra, dec

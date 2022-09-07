# Copyright (C) 2021--2022, the ixpeobssim team.
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

"""Facilities to produce photon lists at the top of the window to be fed into ixpesim.
"""

import numpy
from astropy.io import fits

from ixpeobssim.core.fitsio import xBinTableHDUBase
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.evt.event import xBaseEventList
from ixpeobssim.evt.fmt import COMMON_HEADER_KEYWORDS, xLvl2PrimaryHDU,\
    xBinTableHDUGTI, xBinTableHDURoiTable, xBinTableHDUSpacecraftData,\
    set_telescope_header_keywords, set_time_header_keywords, set_object_header_keywords
from ixpeobssim.instrument.mma import parse_dithering_kwargs
from ixpeobssim.irf.ebounds import ENERGY_GRID
from ixpeobssim.irfgen.du import uv_filter_transparency
from ixpeobssim.irfgen.gpd import window_transparency
from ixpeobssim.irfgen.mma import effective_area, vignetting_spline
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.math_ import modulo_2pi


# pylint: disable=invalid-name, too-many-ancestors, too-many-locals, too-many-arguments



class xBinTableHDUPhotons(xBinTableHDUBase):

    """Binary table description for the PHOTONS extension in ixpesim event lists.
    """

    NAME = 'PHOTONS'
    HEADER_KEYWORDS = COMMON_HEADER_KEYWORDS
    DATA_SPECS = [
        ('SEC'     , 'J', 's'  , 'Integral part of event time (MET)'),
        ('MICROSEC', 'J', 'us' , 'Fractional part of event time (MET)'),
        ('TIME'    , 'D', 's'  , 'Event time in seconds (MET)'),
        ('ENERGY'  , 'E', 'keV', 'Event energy in keV'),
        ('RA'      , 'E', 'deg', 'Right Ascension of the photon'),
        ('DEC'     , 'E', 'deg', 'Declination of the photon'),
        ('DETX'    , 'E', 'mm' , 'X position at the top of the Be window (GPD frame)'),
        ('DETY'    , 'E', 'mm' , 'Y position at the top of the Be window (GPD frame)'),
        ('POL_DEG' , 'E', ''   , 'Polarization degree'),
        ('POL_ANG' , 'E', 'rad', 'Polarization angle'),
        ('SRC_ID'  , 'I', None , 'Monte Carlo source identifier')
    ]



def build_tow_response(irf_set):
    """Build the response at the top of the Be window that is needed to generate the
    photon lists.

    This is assembled "by hand" using the irfgen facilities, and to do this we
    need to know the proper files used to build the response functions in the
    first place. Unfortunately, for historical reasons, these are stored in the
    COMMENT field of the primary header of the FITS files, and therefore we
    have to resort to some string parsing to make this happen. Admittedly, this
    is fragile, and we are guaranteed to break it if we ever change the format of
    the comments. (Note, however, that we are starting with sensible defaults.)
    """
    logger.info('Building the response at the top of the window for DU %d', irf_set.du_id)
    # The effective area we care about is the product of the MMA effectiva area
    # and the transparency of the UV filter.
    aeff = effective_area(ENERGY_GRID, irf_set.du_id)
    aeff *= uv_filter_transparency(ENERGY_GRID)
    # Correct for the different modeling of the Be window assemblies between
    # ixpeobssim and ixpesim.
    scale = window_transparency(ENERGY_GRID) / window_transparency(ENERGY_GRID, None)
    aeff *= scale
    aeff_spline = xInterpolatedUnivariateSpline(ENERGY_GRID, aeff)
    # And we also need the vignetting data.
    vign_spline = vignetting_spline()
    return aeff_spline, vign_spline



class xPhotonList(xBaseEventList):

    """Class describing a photon list.

    This is, in many respects, the equivalent of the evt.event.xEventList
    class, except that it does not represent actual events in the detector, but
    photons at the top of the Be window.
    """

    _TABLE_CLASSES = (xBinTableHDUPhotons, )

    def fill(self, energy, ra, dec, detx, dety, pol_deg, pol_ang):
        """Fill all the relevant columns of the event list.
        """
        self._set_column('ENERGY', energy)
        self._set_column('RA', ra)
        self._set_column('DEC', dec)
        self._set_column('DETX', detx)
        self._set_column('DETY', dety)
        self._set_column('POL_DEG', pol_deg)
        self._set_column('POL_ANG', pol_ang)

    def apply_vignetting(self, vign, ra_pnt, dec_pnt):
        """Apply the effective area vignetting to the event list.
        """
        ra, dec, energy = self.get('RA'), self.get('DEC'), self.get('ENERGY')
        self.apply_vignetting_base(ra, dec, energy, vign, ra_pnt, dec_pnt)

    def apply_fiducial_area(self):
        """Trim out the events falling outside the detector active area.
        """
        detx, dety = self.get('DETX'), self.get('DETY')
        self.apply_fiducial_area_base(detx, dety)

    def _finalize(self, irf_set, **kwargs):
        """Finalize the event list.

        This should be run when there are no more events to be added to the list
        and we're ready to write the list itself to the output file. The method
        basically runs several different task, if and when necessary---in this order:

        * run `apply_fiducial_area()`
        * sort the event list;
        """
        # pylint: disable=unused-argument
        # If there are no events in the list we refrain from doing anything.
        if self.num_events() == 0:
            return
        logger.info('Finalizing the photon list...')
        self.apply_fiducial_area()
        self.sort()
        # Need to rotate the PHI angle by 90 degrees in order to have the origin
        # of the coordinate system for the position angle at the celestial North,
        # see https://bitbucket.org/ixpesw/ixpeobssim/issues/597
        # There has been a lot of back and forth on this one, and we finally
        # convinced ourselves that we need to rotate by 90 degrees here, *and*
        # rotate back by -90 degree in ixpesim for the whole thing to round-trip
        # correctly, i.e., produce the right pattern of Stokes crosstalk, and
        # preserve the instrinsic source pattern in detector coordinates.
        phi = self.get('POL_ANG')
        phi = modulo_2pi(phi + 0.5 * numpy.pi)
        self._set_column('POL_ANG', phi)

    def write_fits(self, creator, roi_model, irf_set, **kwargs):
        """Write the photon list and associated ancillary information to file.
        """
        du_id = irf_set.du_id
        self._finalize(irf_set, **kwargs)
        file_path = kwargs.get('outfile')
        start_met = kwargs.get('start_met')
        duration = kwargs.get('duration')
        stop_met = start_met + duration
        gti_list = kwargs['gti_list']
        ontime = gti_list.total_good_time()

        # Cache the proper args for the common header keywords.
        _time_args = start_met, stop_met, duration, ontime, ontime, 1.
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

        # Create the PHOTONS extension.
        photon_hdu = xBinTableHDUPhotons(self)
        _update_header(photon_hdu)

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
        hdu_list = fits.HDUList([primary_hdu, photon_hdu, gti_hdu, roi_hdu])

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

        # We're good to go.
        hdu_list.info()
        hdu_list.writeto(file_path, overwrite=True)
        hdu_list.close()
        logger.info('Photon list written to %s...', file_path)

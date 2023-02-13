# Copyright (C) 2022, the ixpeobssim team.
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

"""Basic data format definitions for the event lists.
"""

from __future__ import print_function, division

from ixpeobssim.core.fitsio import xPrimaryHDU, xBinTableHDUBase
from ixpeobssim.instrument.du import du_physical_name, du_logical_name
from ixpeobssim.utils.astro import xy_columns_kwargs, set_xy_header_limits, build_wcs
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.time_ import MISSION_START_MJD, MISSION_START_MJDREFF
from ixpeobssim.utils.time_ import met_to_string
from ixpeobssim.utils.units_ import arcsec_to_degrees, degrees_to_arcmin

# pylint: disable=invalid-name, too-many-ancestors, too-many-arguments, no-member


# Header keywords related to the telescope.
_TELESCOPE_HEADER_KEYWORDS = [
    ('TELESCOP', 'IXPE', 'Telescope name'),
    ('INSTRUME', 'GPD', 'Instrument name'),
    ('DETNAM'  , 'N/A', 'Name of the logical detector unit'),
    ('DET_ID'  , 'N/A', 'Name of the physical detector unit')
]

def set_telescope_header_keywords(hdu, du_id):
    """Set the telescope-related header keywords for a given HDU.
    """
    logger.info('Updating telescope header keywords for the %s HDU...', hdu.name)
    header = hdu.header
    header.set('DETNAM', du_logical_name(du_id))
    header.set('DET_ID', du_physical_name(du_id))


# Header keywords related to time.
_TIME_HEADER_KEYWORDS = [
    ('TSTART'  , -1, 'Observation start time in MET'),
    ('TSTOP'   , -1, 'Observation end time in MET'),
    ('DATE-OBS', 'N/A', 'Observation start datetime'),
    ('DATE-END', 'N/A', 'Observation end datetime'),
    ('TELAPSE' , -1, 'TSTOP-TSTART'),
    ('TIMESYS' , 'TT', 'Time system'),
    ('TIMEUNIT', 'seconds', 'Time units'),
    ('TIMEREF' , 'LOCAL', 'Time reference'),
    ('MJDREFI' , MISSION_START_MJD, 'MJD ref day 01 Jan 2017 00:00:00 UTC'),
    ('MJDREFF' , MISSION_START_MJDREFF, 'Frac part of MJD ref (32.184secs+37leapsecs)'),
    ('TIMEZERO', 0., 'Zero time'),
    ('ONTIME'  , -1, 'On source time'),
    ('LIVETIME', -1, 'On source time corrected for dead time'),
    ('DEADC'   , -1, 'Dead time correction'),
    ('DEADAPP' , 'F', 'Has DEADC been applied to data')
]

def set_time_header_keywords(hdu, start_met, stop_met, duration, ontime, livetime, deadc):
    """Set the time-related header keywords for a given HDU.
    """
    logger.info('Updating time header keywords for the %s HDU...', hdu.name)
    header = hdu.header
    header.set('TSTART', start_met)
    header.set('TSTOP', stop_met)
    header.set('DATE-OBS', met_to_string(start_met))
    header.set('DATE-END', met_to_string(stop_met))
    header.set('TELAPSE', duration)
    header.set('ONTIME', ontime)
    header.set('LIVETIME', livetime)
    header.set('DEADC', deadc)


#Header keywords related to the object being observed.
_OBJECT_HEADER_KEYWORDS = [
    ('RA_OBJ' , 'N/A', '[deg] R.A. Object'),
    ('DEC_OBJ', 'N/A', '[deg] Dec Object'),
    ('RA_PNT' , 'N/A', '[deg] RA pointing'),
    ('DEC_PNT', 'N/A', '[deg] Dec pointing'),
    ('OBJECT' , 'N/A', 'Name of observed object')
]

def set_object_header_keywords(hdu, ra_obj, dec_obj, name=None):
    """Set the object-related header keywords for a given HDU.
    """
    logger.info('Updating object header keywords for the %s HDU...', hdu.name)
    header = hdu.header
    header.set('RA_OBJ', ra_obj)
    header.set('DEC_OBJ', dec_obj)
    header.set('RA_PNT', ra_obj)
    header.set('DEC_PNT', dec_obj)
    if name is not None:
        header.set('OBJECT', name)

COMMON_HEADER_KEYWORDS = \
    _TELESCOPE_HEADER_KEYWORDS + _TIME_HEADER_KEYWORDS + _OBJECT_HEADER_KEYWORDS

#Header keywords related to the WCS.
_WCS_HEADER_KEYWORDS = [
    ('EQUINOX' , 2000, 'equinox of celestial coord system'),
    ('RADECSYS', 'FK5', 'celestial coord system')
]

def set_wcs_header_keywords(hdu):
    """Set the WCS-related header keywords for a given HDU.
    """
    logger.info('Updating WCS header keywords for the %s HDU...', hdu.name)
    header = hdu.header
    for key, val, _ in _WCS_HEADER_KEYWORDS:
        logger.debug('Setting %s -> %s', key, val)
        header.set(key, val)


#Header keywords related to file.
_VERSION_HEADER_KEYWORDS = [
    ('FILE_LVL', 'LV2', 'File level'),
    ('LV2_VER' , -1   , 'Version of the LV2 data format')
]

def set_version_keywords(hdu, version):
    """Set the version-related header keywords for a given HDU.
    """
    logger.info('Updating version header keywords for the %s HDU...', hdu.name)
    header = hdu.header
    header.set('LV2_VER', version)


# Basic parameters for the output sky coordinates.
_SKYCOORD_NUM_SIDE_PIXELS = 600
_SKYCOORD_PIXEL_SIZE = arcsec_to_degrees(2.6)


def standard_xy_columns_kwargs(ra0, dec0):
    """Return the appropriate keywords for the X and Y columns in event files.

    See https://bitbucket.org/ixpesw/ixpeobssim/issues/552
    """
    return xy_columns_kwargs(ra0, dec0, _SKYCOORD_NUM_SIDE_PIXELS, _SKYCOORD_PIXEL_SIZE)


def set_standard_xy_header_limits(hdu):
    """Set the TLMIN and TLMAX keywords for the X and Y columns for a given hdu.
    """
    set_xy_header_limits(hdu, _SKYCOORD_NUM_SIDE_PIXELS)


def build_standard_wcs(ra0, dec0):
    """Build the standard WCS object for event files programmatically.

    See https://bitbucket.org/ixpesw/ixpeobssim/issues/552
    """
    return build_wcs(ra0, dec0, _SKYCOORD_NUM_SIDE_PIXELS, _SKYCOORD_PIXEL_SIZE)


# Definition of the WCS origin.
WCS_ORIGIN = 1


def standard_radec_to_xy(ra, dec, ra0, dec0):
    """Convert sky coordinates to X and Y digital coordinates in the sky frame.

    This is builfding on the fly a standard IXPE WCS object centered at the given
    sky position, and using it for the conversion.

    Args
    ----
    ra : array_like
        The array of input Ra coordinates.

    dec : array_like
        The array of input Dec coordinates.

    ra0 : float
        The Ra coordinate of the center of the field in the sky.

    dec0 : float
        The Dec coordinate of the center of the field in the sky.
    """
    return build_standard_wcs(ra0, dec0).wcs_world2pix(ra, dec, WCS_ORIGIN)


def standard_xy_to_radec(x, y, ra0, dec0):
    """Convert X and Y digital coordinates to RA and DEC.

    This is builfding on the fly a standard IXPE WCS object centered at the given
    sky position, and using it for the conversion.

    Args
    ----
    x : array_like
        The array of input X coordinates.

    y : array_like
        The array of input Y coordinates.

    ra0 : float
        The Ra coordinate of the center of the field in the sky.

    dec0 : float
        The Dec coordinate of the center of the field in the sky.
    """
    return build_standard_wcs(ra0, dec0).wcs_pix2world(x, y, WCS_ORIGIN)



class xLvl2PrimaryHDU(xPrimaryHDU):

    """Level 2 primary header definition.
    """

    HEADER_KEYWORDS = [
        ('OBS_ID'  , -1, 'Observation ID'),
        ('CONTNUMB', -1, 'Contact number'),
        ('OBS_MODE', 'OBSERVATION SIMULATION', 'Observation mode'),
        ('SRC_CONF', 'ASTRO', 'Source configuration'),
        ('ORIGIN'  , 'IXPE team', 'Organization responsible for the data'),
        ('CALDBVER', -1, 'CALDB version'),
        ('CLOCKCOR', 'UNKNOWN', 'Whether the time has been corrected')
    ] + COMMON_HEADER_KEYWORDS + _VERSION_HEADER_KEYWORDS



class xBinTableHDUEvents(xBinTableHDUBase):

    """Binary table description for the EVENTS extension of the observation
    output files.
    """

    NAME = 'EVENTS'
    HEADER_KEYWORDS = COMMON_HEADER_KEYWORDS + _WCS_HEADER_KEYWORDS + _VERSION_HEADER_KEYWORDS
    DATA_SPECS = [
        ('TRG_ID'  , 'J', None   , 'Trigger identifier'),
        ('SEC'     , 'J', 's'    , 'Integral part of event time (MET)'),
        ('MICROSEC', 'J', 'us'   , 'Fractional part of event time (MET)'),
        ('TIME'    , 'D', 's'    , 'Event time in seconds (MET)'),
        ('LIVETIME', 'J', 'us'   , 'Live time since the previous event in microseconds'),
        ('PHA'     , 'J', None   , 'Event pulse height'),
        ('PI'      , 'E', None   , 'Event pulse invariant'),
        ('ENERGY'  , 'E', 'keV'  , 'Event energy in keV'),
        ('NUM_CLU' , 'I', ''     , 'Number of clusters in the event'),
        ('DETX'    , 'E', 'mm'   , 'Reconstructed absorption point X (GPD frame)'),
        ('DETY'    , 'E', 'mm'   , 'Reconstructed absorption point Y (GPD frame)'),
        ('RA'      , 'E', 'deg'  , 'Event right ascension'),
        ('DEC'     , 'E', 'deg'  , 'Event declination'),
        ('X'       , 'E', 'pixel', 'Event X position (SKY frame)'),
        ('Y'       , 'E', 'pixel', 'Event Y position (SKY frame)'),
        ('DETPHI'  , 'E', 'rad'  , 'Photolectron emission angle (GPD frame)'),
        ('PHI'     , 'E', 'rad'  , 'Photoelectron emission angle (SKY frame)'),
        ('Q'       , 'E', None   , 'Corrected event q Stokes parameter'),
        ('U'       , 'E', None   , 'Corrected event u Stokes parameter'),
        ('W_MOM'   , 'E', None   , 'Event weight')
    ]

    def __init__(self, ra0, dec0, data=None, keywords=None, comments=None):
        """Overloaded constructor.

        We need this in order to update the DATA_KWARGS class member, so that
        when the constructor of the base class is called, the proper keyword
        arguments are passed to the X and Y FITS column creation.
        """
        xkwargs, ykwargs = standard_xy_columns_kwargs(ra0, dec0)
        self.DATA_KWARGS = {'X': xkwargs, 'Y': ykwargs}
        xBinTableHDUBase.__init__(self, data, keywords, comments)
        # Note that we also have to set the limits for the X and Y columns for the
        # thing to be properly displayed in ds9.
        set_standard_xy_header_limits(self)



class xBinTableHDUMonteCarlo(xBinTableHDUBase):

    """Binary table description for the MONTE_CARLO extension of the observation
    output files, including the additional Monte Carlo fields.
    """

    NAME = 'MONTE_CARLO'
    HEADER_KEYWORDS = [
        ('IRFNAME' , 'N/A', 'name of the IRFs used for the simulation')
    ] + COMMON_HEADER_KEYWORDS
    DATA_SPECS = [
        ('SRC_ID'   , 'I', None     , 'Monte Carlo source identifier'),
        ('MC_ENERGY', 'E', 'keV'    , 'Monte Carlo event energy'),
        ('MC_PHA'   , 'J', None     , 'Monte Carlo pulse height'),
        ('MC_PI'    , 'E', None     , 'Monte Carlo pulse invariant'),
        ('MC_RA'    , 'E', 'degrees', 'Monte Carlo right ascension'),
        ('MC_DEC'   , 'E', 'degrees', 'Monte Carlo declination'),
        ('MC_X'     , 'I', 'degrees', 'Monte Carlo event X position (SKY frame)'),
        ('MC_Y'     , 'I', 'degrees', 'Monte Carlo event Y position (SKY frame)'),
        ('MC_GAIN'  , 'E', None     , 'Relative GEM gain used for the event')
    ]

    def set_irf_name(self, irf_name):
        """Set the IRFNAME keyword.
        """
        self.set_keyword('IRFNAME', irf_name)



class xBinTableHDUGTI(xBinTableHDUBase):

    """Binary table for the good time intervals (GTI).
    """

    NAME = 'GTI'
    HEADER_KEYWORDS = COMMON_HEADER_KEYWORDS
    DATA_SPECS = [
        ('START', 'D', 's', 'GTI start time'),
        ('STOP' , 'D', 's', 'GTI stop time')
    ]


class xBinTableHDURoiTable(xBinTableHDUBase):

    """Binary table for the good time intervals (GTI).
    """

    NAME = 'ROITABLE'
    HEADER_KEYWORDS = [
        ('ROIRA'   , -1, 'right ascension of the ROI center'),
        ('ROIDEC'  , -1, 'declination of the ROI center'),
        ('EQUINOX' , 2000., 'equinox for RA and DEC')
    ] + COMMON_HEADER_KEYWORDS
    DATA_SPECS = [
        ('SRCID'  , 'I'  , None, 'source identifier'),
        ('SRCNAME', 'A20', None, 'source name')
    ]

    def set_center(self, ra0, dec0):
        """Set the keywords for the ROI center.
        """
        self.set_keyword('ROIRA', ra0)
        self.set_keyword('ROIDEC', dec0)



class xBinTableHDUSpacecraftData(xBinTableHDUBase):

    """Binary table for the spacecraft data.
    """

    NAME = 'SC_DATA'
    HEADER_KEYWORDS = [
        ('ROLL'         , -1, 'spacecraft roll angle'),
    ] + COMMON_HEADER_KEYWORDS
    DATA_SPECS = [
        ('MET'          , 'D', 's'      , 'Mission elapsed time'),
        ('RA_PNT'       , 'E', 'degrees', 'Pointing RA'),
        ('DEC_PNT'      , 'E', 'degrees', 'Pointing DEC'),
        ('LAT_GEO'      , 'E', 'degrees', 'Spacecraft latitude'),
        ('LON_GEO'      , 'E', 'degrees', 'Spacecraft longitude'),
        ('ALT_GEO'      , 'E', 'km'     , 'Spacecraft elevation'),
        ('SUN_ANGLE'    , 'E', 'degrees', 'Angle between the Sun and the target'),
        ('IN_SAA'       , 'I', None     , 'SAA flag'),
        ('TARGET_OCCULT', 'I', None     , 'Earth occultation flag')
    ]

    def set_roll_angle(self, roll_angle):
        """Set the roll angle.
        """
        self.set_keyword('ROLL', roll_angle)



class xBinTableHDUTimeline(xBinTableHDUBase):

    """Binary table for the timeline data.
    """

    NAME = 'TIMELINE'
    HEADER_KEYWORDS = COMMON_HEADER_KEYWORDS
    DATA_SPECS = [
        ('START'        , 'D', 's'  , 'Epoch start time'),
        ('STOP'         , 'D', 's'  , 'Epoch stop time'),
        ('IN_SAA'       , 'I', None , 'SAA flag'),
        ('TARGET_OCCULT', 'I', None , 'Earth occultation flag')
    ]



class xBinTableHDUOCTI(xBinTableHDUBase):

    """Binary table for the on-orbit calibration time intervals (OCTIs).
    """

    NAME = 'OCTI'
    HEADER_KEYWORDS = [
        ('CALRUNS', -1, 'Number of on-orbit calibration runs'),
        ('CALTIME', -1, 'Total time for on-orbit calibration in s')
    ] + COMMON_HEADER_KEYWORDS
    DATA_SPECS = [
        ('START', 'D', 's', 'OCTI start time'),
        ('STOP' , 'D', 's', 'OCTI stop time')
    ]

    def set_cal_stats(self, num_runs, total_time):
        """Set the calibration statistics keywords.
        """
        self.set_keyword('CALRUNS', num_runs)
        self.set_keyword('CALTIME', total_time)

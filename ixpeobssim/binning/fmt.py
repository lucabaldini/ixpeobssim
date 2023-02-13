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

"""Format definitions for binned data products.
"""

from __future__ import print_function, division

from ixpeobssim.core.fitsio import xBinTableHDUBase
from ixpeobssim.irf.ebounds import TLMIN, TLMAX


# pylint: disable=invalid-name, too-many-ancestors


class xBinTableHDUPHA1(xBinTableHDUBase):

    """Binary table for binned PHA1 data.
    """

    NAME = 'SPECTRUM'
    HEADER_KEYWORDS = [
        ('HDUCLASS', 'OGIP'),
        ('HDUCLAS1', 'SPECTRUM'),
        ('HDUCLAS2', 'TOTAL'),
        ('HDUCLAS3', 'RATE'),
        ('CHANTYPE', 'PI'),
        ('HDUVERS' , '1.2.1', 'OGIP version number'),
        ('TLMIN1'  , TLMIN, 'first channel number'),
        ('TLMAX1'  , TLMAX, 'last channel number'),
        ('CORRSCAL', 1., 'scaling for correction file'),
        ('POISSERR', False, 'use statistical errors'),
        ('BACKFILE', ''),
        ('CORRFILE', ''),
        ('SYS_ERR' , 0.),
        ('AREASCAL', 1.),
        ('BACKSCAL', 1.)
        ]
    DATA_SPECS = [
        ('CHANNEL' , 'J'),
        ('RATE'    , 'E', 'counts/s'),
        ('STAT_ERR', 'E', 'counts/s'),
        ]



class xBinTableHDUPCUBE(xBinTableHDUBase):

    """Binary table for binned PCUBE data.
    """

    NAME = 'POLARIZATION'
    HEADER_KEYWORDS = []
    # Be careful: if you change any of these, make sure you update the
    # MDP_COLUMNS and POL_COLUMNS class members below, as this all need to be
    # in synch with the MDP and polarization maps and map cubes.
    DATA_SPECS = [
        ('ENERG_LO', 'E', 'keV', 'low energy bound'),
        ('ENERG_HI', 'E', 'keV', 'high energy bound'),
        ('E_MEAN'  , 'E', 'keV', 'average energy within the bin'),
        ('COUNTS'  , 'J', ''   , 'number of counts'),
        ('MU'      , 'E', ''   , 'effective modulation factor'),
        ('W2'      , 'E', ''   , 'sum of weights squared'),
        ('N_EFF'   , 'E', ''   , 'effective number of events w/o acceptance correction'),
        ('FRAC_W'  , 'E', ''   , 'N_EFF / COUNTS'),
        ('MDP_99'  , 'E', ''   , 'minimum detectable polarization at the 99% CL'),
        ('I'       , 'E', ''   , 'I Stokes parameter'),
        ('I_ERR'   , 'E', ''   , '1-sigma uncertainty on I'),
        ('Q'       , 'E', ''   , 'Q Stokes parameter'),
        ('Q_ERR'   , 'E', ''   , '1-sigma uncertainty on Q'),
        ('U'       , 'E', ''   , 'U Stokes parameter'),
        ('U_ERR'   , 'E', ''   , '1-sigma uncertainty on U'),
        ('QN'      , 'E', ''   , 'normalized Q Stokes parameter Q/I'),
        ('QN_ERR'  , 'E', ''   , '1-sigma uncertainty on Q/I'),
        ('UN'      , 'E', ''   , 'normalized U Stokes parameter Q/I'),
        ('UN_ERR'  , 'E', ''   , '1-sigma uncertainty on U/I'),
        ('QUN_COV' , 'E', ''   , 'covariance between QN and UN'),
        ('PD'      , 'E', ''   , 'measured polarization degree'),
        ('PD_ERR'  , 'E', ''   , '1-sigma uncertainty on the polarization degree'),
        ('PA'      , 'E', 'deg', 'measured polarization angle'),
        ('PA_ERR'  , 'E', 'deg', '1-sigma uncertainty on the polarization angle'),
        ('P_VALUE' , 'E', ''   , 'p-value for the null hypothesis (no polarization)'),
        ('CONFID'  , 'E', ''   , 'confidence of the polarization detection'),
        ('SIGNIF'  , 'E', ''   , 'detection significance in equivalent gaussian sigma')
        ]
    COL_NAMES = [col_name for col_name, *_ in DATA_SPECS]
    # Be careful: these need to be in synch with the DATA_SPECS above.
    # MDP_COL_NAMES defines the extensions of the MDP maps and map cubes.
    # POL_COL_NAMES defines the extensions of the polarization maps and map cubes.
    MDP_COL_NAMES = COL_NAMES[2:10]
    POL_COL_NAMES = COL_NAMES[2:]



class xBinTableHDUEBOUNDS(xBinTableHDUBase):

    """Binary table for storing energy bounds.

    .. warning::

       Can we just reuse the same extension from the response matrix?
    """

    NAME = 'EBOUNDS'
    HEADER_KEYWORDS = []
    DATA_SPECS = [
        ('ENERG_LO', 'E', 'keV'),
        ('ENERG_HI', 'E', 'keV')
        ]



class xBinTableHDULC(xBinTableHDUBase):

    """Binary table for binned LC data.
    """

    NAME = 'RATE'
    HEADER_KEYWORDS = []
    DATA_SPECS = [
        ('TIME'    , 'D', 's'     , 'time at the bin center'),
        ('TIMEDEL' , 'D', 's'     , 'bin size'),
        ('EXPOSURE', 'D', 's'     , 'exposure in bin'),
        ('COUNTS'  , 'D', 'counts', 'photon counts'),
        ('ERROR'   , 'E', 'counts', 'statistical errors')
        ]



class xBinTableHDUPP(xBinTableHDUBase):

    """Binary table for binned PP data.
    """

    NAME = 'RATE'
    HEADER_KEYWORDS = []
    DATA_SPECS = [
        ('PHASE'   , 'D', 's'     , 'phase at the bin center'),
        ('PHASEDEL', 'D', 's'     , 'phase bin size'),
        ('COUNTS'  , 'J', 'counts', 'photon counts'),
        ('ERROR'   , 'E', 'counts', 'statistical errors')
        ]



class xBinTableHDUTHETABOUNDS(xBinTableHDUBase):

    """Binary table for storing off-axis angle bounds.
    """

    NAME = 'THETA_BOUNDS'
    DATA_SPECS = [
        ('THETA_LO', 'E', 'arcmin'),
        ('THETA_HI', 'E', 'arcmin')
        ]

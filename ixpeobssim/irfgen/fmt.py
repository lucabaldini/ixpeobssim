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

"""Basic format definitions for the response files.
"""

from __future__ import print_function, division


from ixpeobssim.core.fitsio import xBinTableHDUBase


# pylint: disable=invalid-name, too-many-ancestors


# Specifications for the FITS headers related to the OGIP standards.
OGIP_HEADER_SPECS = [
    ('HDUCLASS', 'OGIP    ', 'format conforms to OGIP standard'),
    ('HDUVERS' , '1.1.0   ', 'Version of format (OGIP memo CAL/GEN/92-002a)'),
    ('HDUDOC'  , 'OGIP memos CAL/GEN/92-002 & 92-002a', 'Documents describing the forma'),
    ('HDUVERS1', '1.0.0   ', 'Obsolete - included for backwards compatibility'),
    ('HDUVERS2', '1.1.0   ', 'Obsolete - included for backwards compatibility')
]


# Additional common header keywords for response files.
RESPONSE_HEADER_SPECS = [
    ('VERSION' , None        , 'Extension version number'),
    ('FILENAME', None        , 'File name'),
    ('CCLS0001', 'CPF'       , 'Dataset is a Calibration Product File'),
    ('CDTP0001', 'DATA'      , 'Calibration file contains data'),
    ('CCNM0001', None        , 'Type of Calibration data'),
    ('CVSD0001', '2017-01-01', 'UTC date when file should first be used'),
    ('CVST0001', '00:00:00'  , 'UTC time when file should first be used'),
    ('CDES0001', None        , 'Description'),
    ('CBD10001', None        , 'Parameter Boundary')
]


class xBinTableHDUMATRIX(xBinTableHDUBase):

    """Binary table for the MATRIX extension of a rmf file.
    """

    NAME = 'MATRIX'
    HEADER_KEYWORDS = [
        ('HDUCLAS1', 'RESPONSE'  , 'dataset relates to spectral response'),
        ('HDUCLAS2', 'RSP_MATRIX', 'dataset is a spectral response matrix'),
        ('CHANTYPE', 'PI '       , 'Detector Channel Type in use (PHA or PI)')
    ] + OGIP_HEADER_SPECS + RESPONSE_HEADER_SPECS
    DATA_SPECS = [
        ('ENERG_LO', 'E', 'keV'),
        ('ENERG_HI', 'E', 'keV'),
        ('N_GRP'   , 'I'),
        ('F_CHAN'  , 'I'),
        ('N_CHAN'  , 'I'),
        ('MATRIX'  , None)
    ]

    def __init__(self, num_channels, data=None, keywords=None, comments=None):
        """Overloaded constructor.

        Note that we have to add the MATRIX column, whose dimension is only
        known at runtime, to the DATA_SPECS
        """
        self.DATA_SPECS[-1] = ('MATRIX', '%dE' % num_channels)
        xBinTableHDUBase.__init__(self, data, keywords, comments)



class xBinTableHDUEBOUNDS(xBinTableHDUBase):

    """Binary table for the MATRIX extension of a rmf file.
    """

    NAME = 'EBOUNDS'
    HEADER_KEYWORDS = [
        ('CHANTYPE', 'PI'             , 'Channel type'),
        ('CONTENT' , 'Response Matrix', 'File content'),
        ('HDUCLAS1', 'RESPONSE'       , 'Extension contains response data  '),
        ('HDUCLAS2', 'EBOUNDS '       , 'Extension contains EBOUNDS')
    ] + OGIP_HEADER_SPECS + RESPONSE_HEADER_SPECS
    DATA_SPECS = [
        ('CHANNEL', 'I'),
        ('E_MIN'  , 'E', 'keV'),
        ('E_MAX'  , 'E', 'keV')
    ]



class xBinTableHDUSPECRESPBase(xBinTableHDUBase):

    """Base class for binary tables describing the SPECRESP extension of a
    response files---depending on the units of the SPECRESP column, that can
    be specifies by overriding the _SPECRESP_UNITS member, this can be
    specialized to effective area/modulation response files or modulation
    factor files.
    """

    NAME = 'SPECRESP'
    HEADER_KEYWORDS = [
        ('HDUCLAS1', 'RESPONSE', 'dataset relates to spectral response'),
        ('HDUCLAS2', 'SPECRESP', 'dataset contains spectral response')
    ] + OGIP_HEADER_SPECS + RESPONSE_HEADER_SPECS
    EBOUND_SPECS = [
        ('ENERG_LO', 'E', 'keV'),
        ('ENERG_HI', 'E', 'keV')
    ]



class xBinTableHDUSPECRESPARF(xBinTableHDUSPECRESPBase):

    """Binary table descriptor for effective area and modulation response files.
    """

    DATA_SPECS = xBinTableHDUSPECRESPBase.EBOUND_SPECS + [
        ('SPECRESP', 'E', 'cm**2')
    ]



class xBinTableHDUSPECRESPMRF(xBinTableHDUSPECRESPARF):

    """Binary table descriptor for effective area and modulation response files.
    """

    pass


class xBinTableHDUSPECRESPMODF(xBinTableHDUSPECRESPBase):

    """Binary table descriptor for modulation factor files.
    """

    DATA_SPECS = xBinTableHDUSPECRESPBase.EBOUND_SPECS + [
        ('SPECRESP', 'E', None)
    ]



class xBinTableHDUVIGNETTING(xBinTableHDUBase):

    """Binary table for the VIGNETTING extension of a arf file.
    """

    NAME = 'VIGNETTING'
    HEADER_KEYWORDS = []
    DATA_SPECS = []

    def __init__(self, data, keywords=None, comments=None):
        """Overloaded constructor.

        Note all the binning is only know at runtime, and the DATA_SPEC,
        accordingly, is entirely assembled within the constructor.
        """
        energy_low, energy_hi, theta, vignetting = data
        num_energy_bins = len(energy_low)
        num_theta_bins = len(theta)
        num_bins = num_energy_bins * num_theta_bins
        assert vignetting.shape == (num_theta_bins, num_energy_bins)
        data = [
            energy_low.reshape((1, num_energy_bins)),
            energy_hi.reshape((1, num_energy_bins)),
            theta.reshape((1, num_theta_bins)),
            vignetting.reshape((1, num_bins))
        ]
        self.DATA_SPECS = [
            ('ENERG_LO'  , '%dE' % num_energy_bins, 'keV'),
            ('ENERG_HI'  , '%dE' % num_energy_bins, 'keV'),
            ('THETA'     , '%dE' % num_theta_bins, 'arcmin'),
            ('VIGNETTING', '%dE' % num_bins)
        ]
        self.HEADER_KEYWORDS = [
            ('TDIM4'     , '(%d, %d)' % (num_energy_bins, num_theta_bins))
        ]
        xBinTableHDUBase.__init__(self, data, keywords, comments)



class xBinTableHDUPSF(xBinTableHDUBase):

    """Binary table for the PSF extension of a psf file.
    """

    NAME = 'PSF'
    DATA_SPECS = [
        ('W'    , 'E', '1/sr'),
        ('SIGMA', 'E', 'arcsec'),
        ('N'    , 'E', '1/sr'),
        ('R_C'  , 'E', 'arcsec'),
        ('ETA'  , 'E')
    ]

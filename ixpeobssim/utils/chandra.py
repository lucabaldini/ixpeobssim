# Copyright (C) 2016--2022 the ixpeobssim team.
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

"""Utilities for the Chandra to IXPE conversion.
"""

from __future__ import print_function, division


import os
import numpy

from astropy.io import fits

from ixpeobssim.utils.logging_ import logger
from ixpeobssim import IXPEOBSSIM_IRF
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline, xInterpolatedBivariateSpline


# pylint: disable=invalid-name, too-many-arguments, no-member


__BASE_IRF_FOLDER = os.path.join(IXPEOBSSIM_IRF, 'fits')


def gti(gtidata):
    """Return the sum of the GTI in the corresponding extension data.
    """
    return (gtidata.field('STOP') - gtidata.field('START')).sum()

def livetime(header):
    """Return the livetime of the observation.
    """
    return header['LIVETIME']

def pointing(header):
    """Return the pointing direction.
    """
    try:
        ra = header['RA_PNT']
        dec = header['DEC_PNT']
    except KeyError:
        ra = header['RA_TARG']
        dec = header['DEC_TARG']
    return ra, dec

def load_arf(detname='ACIS-I'):
    """Load the Chandra effective area data from file.
    """
    detname = detname.lower().replace('-', '_')
    assert detname in ['acis_s', 'acis_i']
    file_path = os.path.join(__BASE_IRF_FOLDER, 'chandra_%s.arf' % detname)
    logger.info('Loading Chandra effective area from %s...', file_path)
    with fits.open(file_path) as hdu_list:
        tbdata = hdu_list['SPECRESP'].data
    _x = 0.5*(tbdata.field('ENERG_LO') + tbdata.field('ENERG_HI'))
    _y = tbdata.field('SPECRESP')
    fmt = dict(xlabel='Energy [keV]', ylabel='Effective area [cm$^2$]')
    return xInterpolatedUnivariateSpline(_x, _y, **fmt)

def load_vign():
    """Load the Chandra vignetting data from file.
    """
    file_path = os.path.join(__BASE_IRF_FOLDER, 'chandra_vignet.fits')
    logger.info('Loading Chandra vignetting from %s...', file_path)
    with fits.open(file_path) as hdu_list:
        tbdata = hdu_list['AXAF_VIGNET'].data
    _x = 0.5*(tbdata.field('ENERG_LO') + tbdata.field('ENERG_HI'))[0,:]
    _y = tbdata.field('THETA')[0,:]
    _vignet = tbdata.field('VIGNET')
    _z = numpy.mean(_vignet[0,:,:,:], axis=0)
    _mask = _y <= 10 #arcmin
    _y = _y[_mask]
    _z = _z[_mask, :]
    fmt = dict(xlabel='Energy [keV]', ylabel='Off-axis angle [arcmin]',
               zlabel='Vignetting')
    return xInterpolatedBivariateSpline(_x, _y, _z.T, **fmt)

def arf_ratio(ixpe_arf, ixpe_vign, chandra_arf, chandra_vign, emin=1., emax=10., thetamax=8.5):
    """Make the ratio between IXPE and Chandra effective area taking in
    account the vignetting.

    This is returning an interpolated bivariate spline in the energy-off-axis
    angle plane (the energy is measured in keV and the off-axis angle in arcmin).
    """
    _x = numpy.arange(emin, emax, 0.005)
    _y = numpy.arange(0., thetamax, 0.05)
    _arf_ratio = ixpe_arf(_x) / chandra_arf(_x)
    _x_mesh, _y_mesh = numpy.meshgrid(_x, _y)
    _vign_ratio = ixpe_vign(_x_mesh, _y_mesh) / chandra_vign(_x_mesh, _y_mesh)
    _z = _vign_ratio * _arf_ratio
    fmt = dict(xlabel='Energy [keV]', ylabel='Off-axis angle [arcmin]',
        zlabel='Effective area ratio')
    return xInterpolatedBivariateSpline(_x, _y, _z.T, **fmt)

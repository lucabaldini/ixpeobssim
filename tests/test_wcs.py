#!/usr/bin/env python
#
# Copyright (C) 2019, the ixpeobssim team.
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

from __future__ import print_function, division

import os
import unittest
import sys

import numpy
from astropy import wcs

from ixpeobssim.srcmodel.polarization import xStokesSkyMap
from ixpeobssim.binning.base import xEventBinningBase
from ixpeobssim.core.fitsio import xFITSImageBase
from ixpeobssim import IXPEOBSSIM_CONFIG_FITS
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()


class testWCS(unittest.TestCase):

    """Unit test for handling WCS.
    """

    @staticmethod
    def _file_path(file_name):
        """
        """
        return os.path.join(IXPEOBSSIM_CONFIG_FITS, file_name)

    def test(self):
        """
        """
        ra = 350.866590103
        dec = 58.8117384439
        delta = 0.00279893338726
        nside = 50
        wcs_ = wcs.WCS(naxis=2)
        wcs_.wcs.crpix = [25., 25.]
        wcs_.wcs.cdelt = numpy.array([-delta, delta])
        wcs_.wcs.crval = [ra, dec]
        wcs_.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        wcs_.pixel_shape = (nside, nside)
        cx, cy = xStokesSkyMap.wcs_center(wcs_)
        self.assertAlmostEqual(cx, ra)
        self.assertAlmostEqual(cy, dec)

    def test_msh1552(self):
        """Test on the MSH 1552 models provided by Niccolo` Bucciantini.

        By visual inspection on ds9 the center of the map is approximately at
        (228.481, -59.135), and the half width of the field is 0.1 degrees.
        """
        ra, dec = 228.481, -59.135
        pd_file_path = self._file_path('MSH_Rad_PF.fits')
        pa_file_path = self._file_path('MSH_Rad_PA.fits')
        radius = 0.1
        map_ = xStokesSkyMap.load_from_pda(pd_file_path, pa_file_path)
        cx, cy = map_.center
        self.assertAlmostEqual(cx, ra, places=2)
        self.assertAlmostEqual(cy, dec, places=2)
        print('MSH1552 radius: ', map_.radius)
        #self.assertAlmostEqual(map_.radius, radius, places=2)

    @unittest.skip('Initialization from xy deprecated')
    def test_tycho(self):
        """
        By visual inspection on ds9 the center of the map is approximately at
        (6.33, 64.14), and the half width of the field is 0.18 degrees
        """
        ra, dec = 6.33, 64.14
        radius = 0.18
        px_file_path = self._file_path('polx_0.4_pf_0.90_radial.fits')
        py_file_path = self._file_path('poly_0.4_pf_0.90_radial.fits')
        map_ = xStokesSkyMap.load_from_xy(px_file_path, py_file_path)
        cx, cy = map_.center
        self.assertAlmostEqual(cx, ra, places=2)
        self.assertAlmostEqual(cy, dec, places=2)
        print('Tycho radius: ', map_.radius)
        #self.assertAlmostEqual(map_.radius, radius, places=2)

    def test_binning_wcs(self, ra0=30., dec0=45., npix=100):
        """This is meant to debug the binning WCS-related routines.
        """
        kwargs = dict(xref=ra0, yref=dec0, npix=npix, proj='TAN')
        wcs_ = xEventBinningBase._build_image_wcs(**kwargs)
        print(wcs_)
        size = 100000
        sigma = 0.0075
        ra = numpy.random.normal(ra0, sigma / numpy.cos(numpy.radians(dec0)), size)
        dec = numpy.random.normal(dec0, sigma, size)
        x, y, binning = xEventBinningBase._pixelize_skycoords(ra, dec, wcs_)
        # Confirm that any offset in the open interval (0, 1) gives the same
        # answer.
        for offset in [0.01, 0.99]:
            _x, _y, _ = xEventBinningBase._pixelize_skycoords(ra, dec, wcs_, offset)
            self.assertTrue(numpy.allclose(_x, x))
            self.assertTrue(numpy.allclose(_y, y))
        data, _, _ = numpy.histogram2d(x, y, bins=binning)
        xFITSImageBase.make_plot(data, wcs_)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

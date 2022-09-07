#!/usr/bin/env python
#
# Copyright (C) 2020, the ixpeobssim team.
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

"""Unit tests for the model intensity facilities.
"""

from __future__ import print_function, division

import unittest
import sys
import os

from ixpeobssim.srcmodel.roi import xPointSource, xUniformDisk, xGaussianDisk
from ixpeobssim.srcmodel.roi import xUniformAnnulus, xExtendedSource
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.units_ import arcmin_to_degrees, degrees_to_arcmin
from ixpeobssim.core.fitsio import xFITSImageBase
from ixpeobssim.utils.astro import build_wcs
from ixpeobssim import IXPEOBSSIM_CONFIG
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()


# pylint: disable=invalid-name


class TestModelIntensity(unittest.TestCase):

    """Unit test for model broadband polarization facilities.
    """

    @classmethod
    def setUpClass(cls):
        """Setup the test.
        """
        spec = power_law(1., 2.)
        pol_deg = constant(0.)
        pol_ang = constant(0.)
        cls.source_args = spec, pol_deg, pol_ang

    def test_intensity_base(self, side=12., nside=101):
        """Unit test for all the source classes providing an analytical
        implementation of build_intensity_map().
        """
        ra, dec = 45., 45.
        r1 = arcmin_to_degrees(2.5)
        r2 = arcmin_to_degrees(5.)
        # We offset the wcs to highlight possible issue with data transposition.
        pixel_size = arcmin_to_degrees(side) / nside
        wcs_ = build_wcs(ra, dec + 5 * pixel_size, nside, pixel_size)
        point_source = xPointSource('Point source', ra, dec, *self.source_args)
        uniform_disk = xUniformDisk('Uniform disk', ra, dec, r1, *self.source_args)
        gauss_disk = xGaussianDisk('Gaussian disk', ra, dec, r1, *self.source_args)
        annulus = xUniformAnnulus('Annulus', ra, dec, r1, r2, *self.source_args)
        for source in [point_source, uniform_disk, gauss_disk, annulus]:
            data = source.build_intensity_map(wcs_)
            self.assertAlmostEqual(data.sum(), 1.)
            plt.figure('%s intensity map' % source.name)
            xFITSImageBase.make_plot(data, wcs_, zlabel='Intensity [a. u.]')

    def test_extended(self, nside=101):
        """Test for extended sources.
        """
        file_path = os.path.join(IXPEOBSSIM_CONFIG, 'fits', 'casa_1p5_3p0_keV.fits')
        source = xExtendedSource('Cas A', file_path, *self.source_args)
        ra, dec = source.image.center()
        pixel_size = 2. * source.image.inner_radius() / nside
        #side = 2. * degrees_to_arcmin(source.image.inner_radius())
        wcs_ = build_wcs(ra, dec, nside, pixel_size)
        data = source.build_intensity_map(wcs_)
        self.assertAlmostEqual(data.sum(), 1.)
        plt.figure('%s intensity map' % source.name)
        xFITSImageBase.make_plot(data, wcs_, zlabel='Intensity [a. u.]')
        plt.figure('%s original template' % source.name)
        source.image.plot()



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

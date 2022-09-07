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

"""Unit tests for the misc.astro module.
"""


from __future__ import print_function, division

import unittest

import numpy
from astropy.coordinates import SkyCoord

from ixpeobssim.irf import load_psf
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.astro import angular_separation, build_wcs, wcs_digitize
from ixpeobssim.utils.units_ import arcmin_to_degrees


# pylint: disable=invalid-name


class TestAstro(unittest.TestCase):

    """Unit test for misc.astro module.
    """

    def test_angular_separation(self):
        """Test the angular separation.
        """
        ra0, dec0 = 30., 45.
        ref = SkyCoord(ra0, dec0, unit='deg')
        psf = load_psf()
        ra, dec = psf.smear_single(ra0, dec0, 10000)
        skycoord_sep = SkyCoord(ra, dec, unit='deg').separation(ref).degree
        astro_sep = angular_separation(ra, dec, ra0, dec0)
        self.assertTrue(numpy.allclose(skycoord_sep, astro_sep))

    def test_wcs(self, ra0=30., dec0=45., wcs_side=15., nside=3):
        """Test the WCS facilities.
        """
        pixel_size = arcmin_to_degrees(wcs_side) / nside
        wcs_ = build_wcs(ra0, dec0, nside, pixel_size)
        logger.info(wcs_)
        x, y = wcs_.wcs_world2pix(ra0, dec0, 1)
        self.assertAlmostEqual(float(x), (nside + 1) / 2.)
        self.assertAlmostEqual(float(y), (nside + 1) / 2.)
        x, y = wcs_.wcs_world2pix(ra0, dec0, 0)
        ra, dec = tuple([float(val) for val in wcs_.wcs_pix2world(x, y, 0)])
        self.assertAlmostEqual(ra, ra0)
        self.assertAlmostEqual(dec, dec0)
        n = 10
        target = numpy.zeros((nside, nside))
        idx = int((nside - 1) / 2)
        target[idx, idx] = float(n)
        data = wcs_digitize(wcs_, numpy.full(n, ra0), numpy.full(n, dec0))
        logger.info(data)
        self.assertTrue(numpy.allclose(data, target))



if __name__ == '__main__':
    unittest.main()

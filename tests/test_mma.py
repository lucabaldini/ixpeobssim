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

from __future__ import print_function, division

import unittest
import numpy

from ixpeobssim.instrument.mma import _sky_to_gpd_naive, _gpd_to_sky_naive
from ixpeobssim.instrument.mma import sky_to_gpd, gpd_to_sky
from ixpeobssim.instrument import DU_IDS


class TestMma(unittest.TestCase):

    """Unit test for the instrument.mma module.
    """

    @classmethod
    def setUpClass(cls):
        """
        """
        cls.ra0 = 30.
        cls.dec0 = 30.
        grid = numpy.linspace(-0.05, 0.05, 11)
        cls.ra = cls.ra0 + grid
        cls.dec = cls.dec0 + grid

    def test_naive(self):
        """Test the "naive" transformation.
        """
        # The center of the ROI, with no offset, must map into the origin of the
        # focal plane
        x, y = _sky_to_gpd_naive(self.ra0, self.dec0, self.ra0, self.dec0)
        self.assertAlmostEqual(x, 0.)
        self.assertAlmostEqual(y, 0.)
        # Test the transformations are reversible.
        x, y = _sky_to_gpd_naive(self.ra, self.dec, self.ra0, self.dec0)
        ra, dec = _gpd_to_sky_naive(x, y, self.ra0, self.dec0)
        self.assertTrue(numpy.allclose(self.ra, ra))
        self.assertTrue(numpy.allclose(self.dec, dec))

    def test_no_dither(self):
        """Test the actual transformations without dithering.
        """
        roll_angle = 0.
        for du_id in DU_IDS:
            detx, dety = sky_to_gpd(self.ra, self.dec, None, self.ra0, self.dec0,
                                    du_id, roll_angle)
            ra, dec = gpd_to_sky(detx, dety, None, self.ra0, self.dec0,
                                 du_id, roll_angle)
            self.assertTrue(numpy.allclose(self.ra, ra))
            self.assertTrue(numpy.allclose(self.dec, dec))

    def test_complete(self):
        """Test the full damned thing.
        """
        roll_angle = 10.
        dither_params = (1.6, 907., 101., 449.)
        t = numpy.random.uniform(0., 100000., size=len(self.ra))
        for du_id in DU_IDS:
            detx, dety = sky_to_gpd(self.ra, self.dec, t, self.ra0, self.dec0,
                                    du_id, roll_angle, dither_params)
            ra, dec = gpd_to_sky(detx, dety, t, self.ra0, self.dec0,
                                 du_id, roll_angle, dither_params)
            self.assertTrue(numpy.allclose(self.ra, ra))
            self.assertTrue(numpy.allclose(self.dec, dec))



if __name__ == '__main__':
    unittest.main()

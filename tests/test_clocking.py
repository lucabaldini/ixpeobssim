#!/usr/bin/env python
#
# Copyright (C) 2018, the ixpeobssim team.
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

from ixpeobssim.instrument.gpd import phi_to_detphi, detphi_to_phi, rotate_detxy
from ixpeobssim.instrument.du import du_rotation_angle
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt


class TestClocking(unittest.TestCase):

    """Unit test for the DU clocking
    """

    def test_xy_rotation(self):
        """
        """
        x = 1.
        y = 0.
        for du_id in DU_IDS:
            rx, ry = rotate_detxy(x, y, du_id)
            logger.info('(1, 0) on focal plane --> (%f, %f) on DU%i' %\
                (rx, ry, du_id))
            rotation_angle = numpy.arctan2(ry, rx)
            self.assertAlmostEqual(du_rotation_angle(du_id), rotation_angle)


    def test_phi_rotation(self, roll_angle=0.):
        """
        """
        phi = 0.
        for du_id in DU_IDS:
            rotation_angle = phi_to_detphi(phi, du_id, roll_angle)
            logger.info('%.1f deg on focal plane --> %.1f deg on DU%i: ' %\
                (phi, numpy.degrees(rotation_angle), du_id))
            self.assertAlmostEqual(du_rotation_angle(du_id), rotation_angle)
            phi2 = detphi_to_phi(rotation_angle, du_id, roll_angle)
            self.assertAlmostEqual(phi2, phi)



if __name__ == '__main__':
    unittest.main()

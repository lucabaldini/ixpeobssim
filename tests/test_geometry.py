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
import sys

import numpy

from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.core.geometry import xPoint, xLine, xRay


if sys.flags.interactive:
    plt.ion()


class TestGeometry(unittest.TestCase):

    """Unit test for the core.geometry module.
    """

    def test(self):
        """Test the native O table.
        """
        p = xPoint(0., 0., 0)
        self.assertAlmostEqual(p.norm(), 0.)
        print(p)
        r = xRay(p, 45., 45.)
        print(r)
        print(r.xdir, r.ydir, r.zdir)
        norm = r.xdir**2. + r.ydir**2. + r.zdir**2.
        self.assertAlmostEqual(norm, 1.)
        print(r.xintersect(10.))
        print(r.yintersect(10.))
        print(r.zintersect(10.))
        r = xRay(p, 0., 0.)
        print(r)
        print(r.xdir, r.ydir, r.zdir)
        norm = r.xdir**2. + r.ydir**2. + r.zdir**2.
        self.assertAlmostEqual(norm, 1.)
        print(r.xintersect(0.))
        print(r.yintersect(0.))
        print(r.zintersect(0.))
        print(r.xintersect(10.))
        print(r.yintersect(10.))
        print(r.zintersect(10.))
        self.assertFalse(p.unphysical())
        self.assertTrue(xPoint.unphysical_point().unphysical())
        p1 = xPoint(1., 0., 1.)
        p2 = xPoint(0., 1., 1.)
        print(p1 + p2)
        print(p1 - p2)
        self.assertFalse(p1 == p2)
        l = xLine(p1, p2)
        print(l, l.length())



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

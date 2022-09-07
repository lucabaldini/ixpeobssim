#!/usr/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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


"""Unit test for piecewise splines.
"""


import sys

import unittest
import numpy

from ixpeobssim.core.spline import xInterpolatedPiecewiseUnivariateSpline
from ixpeobssim.utils.matplotlib_ import plt, setup_gca

if sys.flags.interactive:
    plt.ion()


class TestPiecewiseSpline(unittest.TestCase):

    """Unit test for univariate splines.
    """

    def test_base(self):
        """
        """
        x = numpy.linspace(1., 12., 25)
        y = numpy.sqrt(x) + 2. * (x > 4.5) + 3. * (x > 8.5)
        fmt = dict(xlabel='x [a. u.]', ylabel='y [a. u.]')
        spline = xInterpolatedPiecewiseUnivariateSpline(x, y, (4.5, 8.5), **fmt)
        for _x, _y in zip(x, y):
            self.assertAlmostEqual(spline(_x), _y)
        self.assertTrue(numpy.allclose(spline(x), y))
        plt.figure('Piecewise spline')
        spline.plot(overlay=True)
        setup_gca(grids=True)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

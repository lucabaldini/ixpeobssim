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

from __future__ import print_function, division


"""Unit test for the core.spline module.
"""

import sys

import unittest
import numpy

from ixpeobssim.core.spline import xStepFunction

from ixpeobssim.utils.matplotlib_ import plt, last_line_color

if sys.flags.interactive:
    plt.ion()


class TestStepFunction(unittest.TestCase):

    """Unit test for the step function.
    """

    def test(self):
        """
        """
        x = numpy.linspace(0., 10., 11)
        y = x[:-1]
        f = xStepFunction(x, y, xlabel='x [a. u.]', ylabel='y [a. u.]')
        plt.figure('Simple step function')
        f.plot(annotate=True)
        x = numpy.linspace(0., 10., 20)
        plt.plot(x, f(x), 'o')



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

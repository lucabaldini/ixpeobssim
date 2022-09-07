#!/usr/bin/env python
#
# Copyright (C) 2016--2018, the ixpeobssim team.
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


"""Unit tests for the irfgen.constants module.
"""


import unittest
import sys

import numpy

from ixpeobssim.utils.units_ import atm_to_mbar
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.irfgen.constants import *

if sys.flags.interactive:
    plt.ion()


class testConstants(unittest.TestCase):

    """Unit test for the irfgen.constants module.
    """

    def test_dme_density(self):
        """
        """
        self.assertAlmostEqual(dme_density(0., atm_to_mbar(1.)), DME_REF_DENSITY)
        self.assertAlmostEqual(dme_density_perfect(0., atm_to_mbar(1.)), DME_REF_DENSITY)
        d1 = dme_density_perfect(20., 800.)
        d2 = dme_density(20., 800.)
        delta = (d1 - d2) / d2
        # This is from the dme test in gpdsw.
        self.assertAlmostEqual(delta, 0.00711, 2)
        plt.figure('DME density scaling')
        T = numpy.linspace(0., 30., 100)
        p = 800.
        plt.plot(T, 1000. * dme_density_perfect(T, p), label='Perfect gas scaling')
        plt.plot(T, 1000. * dme_density(T, p), label='Custom scaling')
        setup_gca(xlabel='Temperature [$^\\circ$C]', ylabel='Density mg cm$^{-3}$',
                  legend=True, grids=True)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

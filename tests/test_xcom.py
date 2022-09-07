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

import unittest
import sys

import numpy

from ixpeobssim.irfgen.xcom import xCrossSectionTable
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()



class TestXcom(unittest.TestCase):

    """Unit test for the XCOM interface.
    """
    
    def test_plot(self):
        """
        """
        for identifier in ['Al', 'Be', 'Cu', 'He']:
            table = xCrossSectionTable(identifier)
            table.plot()

    def test_dme_efficiency(self, thickness=1.0, density=1.69e-3):
        """
        """
        plt.figure('DME photoabsorption efficiency')
        table = xCrossSectionTable('DME')
        energy = numpy.linspace(1, 15, 100)
        eff = table.photoabsorption_efficiency(energy, thickness, density)
        fmt = dict(xlabel='Energy [keV]', ylabel='Absorption efficiency')
        spline = xInterpolatedUnivariateSplineLinear(energy, eff, **fmt)
        spline.plot(logx=True, logy=True)
        plt.axis([energy.min(), energy.max(), None, None])



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

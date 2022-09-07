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
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.irfgen.astar import xAlphaStoppingPowerTable
from ixpeobssim.irfgen.astar import load_alpha_stopping_power_data
from ixpeobssim.irfgen.astar import load_dme_alpha_stopping_power_data
from ixpeobssim.irfgen.constants import BE_DENSITY, AL_DENSITY

if sys.flags.interactive:
    plt.ion()



class TestAstar(unittest.TestCase):

    """Unit test for the astar module in irfgen.
    """

    def test_be(self):
        """Test the native Be table.
        """
        table = load_alpha_stopping_power_data('Be', BE_DENSITY)
        plt.figure('Be stopping power')
        table.stopping_power.plot(logx=True, logy=True, grids=True)
        plt.figure('Be range')
        table.csda_range.plot(logx=True, logy=True, grids=True)
        setup_gca(xmin=0.01, xmax=10., ymin=1.e-5, ymax=1.e-2)
        plt.figure('Be Bragg curve')
        table.bragg_curve().plot(grids=True)
        energy = 10.
        energy_profile = table._energy_profile(energy)
        plt.figure('Be energy profile @ %.2f MeV' % energy)
        energy_profile.plot(grids=True)

    def _test_al(self):
        """Test the native Al table.
        """
        table = load_alpha_stopping_power_data('Al', AL_DENSITY)
        plt.figure('Al stopping power')
        table.stopping_power.plot(logx=True, logy=True, grids=True)
        plt.figure('Al range')
        table.csda_range.plot(logx=True, logy=True, grids=True)
        setup_gca(xmin=0.01, xmax=10., ymin=1.e-5, ymax=1.e-2)
        plt.figure('Al Bragg curve')
        table.bragg_curve().plot(grids=True)

    def test_dme(self, temperature=20., pressure=0.789539):
        """Test the actual DME calculation.
        """
        table = load_dme_alpha_stopping_power_data()
        plt.figure('DME stopping power')
        table.stopping_power.plot(logx=True, logy=True, grids=True)
        plt.figure('DME range')
        table.csda_range.plot(logx=True, logy=True, grids=True)
        setup_gca(xmin=0.01, xmax=10., ymin=1.e-2, ymax=10.)
        plt.figure('DME Bragg curve')
        table.bragg_curve().plot(grids=True)
        energy = 10.
        energy_profile = table._energy_profile(energy)
        plt.figure('DME energy profile @ %.2f MeV' % energy)
        energy_profile.plot(grids=True)
        self.assertAlmostEqual(energy_profile(0.), energy)
        self.assertAlmostEqual(energy_profile(table.csda_range(energy)), 0.)
        self.assertAlmostEqual(energy_profile(100.), 0.)
        energy_loss = table.energy_loss(energy)
        self.assertAlmostEqual(energy_loss(10., 0.), 0.)
        self.assertAlmostEqual(energy_loss(10., 1000.), 10.)
        print(energy_loss(10., 1.))
        print(energy_loss(5., 1.))
        print(energy_loss(3., 1.))
        print(energy_loss(2., 1.))
        print(energy_loss(1., 1.))


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

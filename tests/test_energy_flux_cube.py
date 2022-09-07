#!/usr/bin/env python
#
# Copyright (C) 2016, the ixpeobssim team.
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

from ixpeobssim.instrument.charging import charging_tau, effective_tau
from ixpeobssim.instrument.charging import asymptotic_gain_delta, calculate_gain
from ixpeobssim.instrument.charging import gain_profile, xEnergyFluxCube
from ixpeobssim.utils.matplotlib_ import plt, last_line_color

if sys.flags.interactive:
    plt.ion()


class TestEnergyFluxCube(unittest.TestCase):

    """Unit test for the GEM charging facilities.
    """

    def test_uniform(self, tmax=20000., num_events=1000000, half_size=1., energy=2.7):
        """
        """
        t = numpy.linspace(0, tmax, num_events)
        x = numpy.random.uniform(-half_size, half_size, size=num_events)
        y = numpy.random.uniform(-half_size, half_size, size=num_events)
        E = numpy.full(t.shape, energy)
        tbinning = numpy.linspace(0, tmax, 500)
        cube = xEnergyFluxCube(200, tbinning)
        cube.fill(x, y, t, E)

        plt.figure('Average energy flux')
        cube.plot()

        plt.figure('Gain')
        t = numpy.linspace(0, tmax * 0.99, 250)
        for _x, _y in ((0., 0.), (3., 3.)):
            x = numpy.full(t.shape, _x)
            y = numpy.full(t.shape, _y)
            plt.plot(t, cube.gain(x, y, t), label='(%.3f, %.3f)' % (_x, _y))
        energy_flux = num_events / tmax * energy / 4 * half_size**2.
        profile = gain_profile(energy_flux)
        plt.plot(t, profile(t), ls='dashed', label='Analytical calculation')
        plt.axis([None, None, 0.92, 1.01])
        plt.legend()
        plt.grid(which='both')
        plt.xlabel('Time [s]')
        plt.ylabel('Ralative gain')

        plt.figure('Gain map')
        cube.gain_map(10).plot()

        plt.figure('Time slice')
        cube.time_slice(10, 100).plot()



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

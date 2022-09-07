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
from ixpeobssim.utils.matplotlib_ import plt, last_line_color

if sys.flags.interactive:
    plt.ion()


class TestCharging(unittest.TestCase):

    """Unit test for the GEM charging facilities.
    """

    def test_model_parameters(self):
        """
        """
        source_labels = [('10 mCrab point source', 0.19),
                         ('Cas A', 0.42),
                         ('100 mCrab point source', 1.9),
                         ('Crab', 13.)
                         ]

        r = numpy.logspace(-2, 2)

        plt.figure('Charging time constant')
        plt.plot(r, effective_tau(r))
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(which='both')
        plt.axis([None, None, 1e3, 2e5])
        plt.xlabel('Energy flux [kev mm$^{-2}$ s$^{-1}$]')
        plt.ylabel('Effective $\\tau$ [s]')
        for l, x in source_labels:
            y = effective_tau(x)
            plt.plot(x, y, 'o', color=last_line_color())
            fmt = dict(rotation=15., ha='left', va='bottom')
            plt.text(x, y, l, color=last_line_color(), **fmt)

        plt.figure('Asymptotic gain delta')
        plt.plot(r, asymptotic_gain_delta(r))
        plt.xscale('log')
        plt.grid(which='both')
        plt.xlabel('Energy [kev mm$^{-2}$ s$^{-1}$]')
        plt.ylabel('Asymptotic gain variation $\\delta$')
        for l, x in source_labels:
            y = asymptotic_gain_delta(x)
            plt.plot(x, y, 'o', color=last_line_color())
            fmt = dict(rotation=15., ha='left', va='bottom')
            plt.text(x, y, l, color=last_line_color(), **fmt)

    def test_gain_calculation(self):
        """
        """
        time_ = numpy.linspace(0., 5e5, 10000)
        plt.figure('Simulated gain profile')
        for r in [0.1, 1., 10.]:
            energy_flux = numpy.full(time_.shape, r)
            gain = calculate_gain(time_, energy_flux)
            plt.plot(time_, gain, label='%.1f keV mm$^{-2} s^{-1}$' % r)
            #print(gain[-1] - 1., asymptotic_gain_delta(r))
        plt.grid(which='both')
        plt.axis([None, None, 0.9, 1.025])
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Relative gain')

    def test_calibration(self):
        """Setup a plausible bench to test the effect.
        """
        area = 25. # mm2
        step = 10. # s
        # Measurement setup: list of (energy, rate, duration) tuples.
        setup = [(2.7, 150., 100000.), (2.7, 10, 500000.), (2.7, 0.1, 800000.),
                 (6.4, 150., 50000.), (6.4, 10, 250000.), (6.4, 0.1, 890000.)]
        tstart = 0.
        time_ = numpy.array([])
        energy_flux = numpy.array([])
        for energy, rate, duration in setup:
            tstop = tstart + duration
            t = numpy.arange(tstart, tstop, step)
            r = energy * rate / area
            time_ = numpy.append(time_, t)
            energy_flux = numpy.append(energy_flux, numpy.full(t.shape, r))
            tstart = tstop
        gain = calculate_gain(time_, energy_flux)
        time_ /= 86400.

        plt.figure('Input energy flux')
        plt.plot(time_, energy_flux)
        plt.grid(which='both')
        plt.xlabel('Time [days]')
        plt.ylabel('Energy flux [kev mm$^{-2}$ s$^{-1}$]')
        plt.yscale('log')

        plt.figure('Gain')
        plt.plot(time_, gain)
        plt.grid(which='both')
        plt.xlabel('Time [days]')
        plt.ylabel('Relative gain')
        plt.axis([None, None, None, 1.01])

        t = 0.
        for energy, rate, duration in setup:
            g = gain[int(t / step)]
            plt.plot(t / 86400., g, 'o', color='orange')
            plt.text(t / 86400., g , ' %.1f keV @ %.1f Hz (on %d mm$^2$)' %\
                    (energy, rate, area), color='orange', size='small')
            t += duration



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

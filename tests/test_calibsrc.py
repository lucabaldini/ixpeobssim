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

import unittest
import sys

import numpy

from ixpeobssim.core.hist import xHistogram1d, xGpdMap2d
from ixpeobssim.irf import load_irf_set
from ixpeobssim.srcmodel.calibsrc import xCalC, xMonochromaticUnpolarizedFlatField
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca

if sys.flags.interactive:
    plt.ion()


class TestCalibrationSources(unittest.TestCase):

    """Unit test for the calibration sources.
    """

    @classmethod
    def setUpClass(cls):
        """Set up the test.
        """
        cls.irf_set = load_irf_set()

    def test_cal_C(self, rate=100., duration=10000., deadtime=0.00108):
        """Test Cal C onboard sources.
        """
        calC = xCalC(rate)
        print(calC)
        kwargs = dict(start_met=0., duration=duration, deadtime=deadtime)
        event_list = calC.rvs_event_list(self.irf_set, **kwargs)
        # Diagnostic plots.
        plt.figure('Cal C energy spectrum')
        energy = event_list.energy()
        binning = numpy.linspace(0., 10., 200)
        hist = xHistogram1d(binning, xlabel='Energy [keV]').fill(energy)
        hist.plot()
        plt.figure('Cal C event times')
        dt = numpy.diff(event_list.time())
        binning = numpy.linspace(0., 10. / rate , 200)
        hist = xHistogram1d(binning, xlabel='Delta event time [s]').fill(dt)
        plt.yscale('log')
        hist.plot()
        plt.figure('Cal C morphology')
        x, y = event_list.detector_coordinates()
        hist = xGpdMap2d(100).fill(x, y)
        hist.plot()

    def test_flat_field(self, rate=100., duration=10000., deadtime=0.00108):
        """Test a flat field.
        """
        flat_field = xMonochromaticUnpolarizedFlatField(rate, 2.7)
        print(flat_field)
        kwargs = dict(start_met=0., duration=duration, deadtime=deadtime)
        event_list = flat_field.rvs_event_list(self.irf_set, **kwargs)
        plt.figure('Flat-field energy spectrum')
        energy = event_list.energy()
        binning = numpy.linspace(0., 10., 200)
        hist = xHistogram1d(binning, xlabel='Energy [keV]').fill(energy)
        hist.plot()
        plt.figure('Flat-field morphology')
        x, y = event_list.detector_coordinates()
        hist = xGpdMap2d(100).fill(x, y)
        hist.plot()



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

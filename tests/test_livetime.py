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


"""Unit test for the LIVETIME column added in conjunction with issue #392.
"""

import unittest
import sys
import os

import numpy
from astropy.io import fits

from ixpeobssim.evt.event import xEventFile
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt
import ixpeobssim.core.pipeline as pipeline

if sys.flags.interactive:
    plt.ion()


class TestLivetime(unittest.TestCase):

    """Unit test for the event list livetime.
    """

    def test_livetime(self, deadtime=0.00108):
        """
        """
        pipeline.reset('toy_point_source', deadtime=deadtime, overwrite=True)
        for file_path in pipeline.xpobssim(duration=10000):
            file_name = os.path.basename(file_path)
            event_file = xEventFile(file_path)
            met = event_file.time_data()
            livetime = event_file.livetime_data() * 1.e-6
            start_met = event_file.min_good_time()
            stop_met = event_file.max_good_time()
            # The livetime for the first event must be equal (up to the microsecond)
            # to the time difference between the first event and the start of
            # the observation.
            self.assertAlmostEqual(livetime[0], met[0] - start_met, places=5)
            self.assertAlmostEqual(livetime.sum(), event_file.livetime())
            self.assertTrue((livetime >= 0.).all())
            self.assertTrue(numpy.allclose(numpy.diff(met) - deadtime, livetime[1:], atol=1.e-5))
            plt.figure('Delta event time')
            dt = numpy.diff(met)
            binning = numpy.linspace(0., dt.max(), 100)
            hist = xHistogram1d(binning).fill(dt)
            hist.plot(label=file_name)
            plt.figure('Livetime')
            binning = numpy.linspace(0., livetime.max(), 100)
            hist = xHistogram1d(binning).fill(livetime)
            hist.plot(label=file_name)
        plt.figure('Delta event time')
        plt.yscale('log')
        plt.legend()
        plt.figure('Livetime')
        plt.yscale('log')
        plt.legend()



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

#!/urs/bin/env python
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

from ixpeobssim.instrument.fcw import xOnOrbitCalibrationRun, xOnOrbitCalibrationPattern



class TestFcw(unittest.TestCase):


    def test_run(self, start_met=0., stop_met=1000.):
        """
        """
        run = xOnOrbitCalibrationRun(start_met, stop_met, 1, None)
        print(run)

    def test_pattern(self):
        """
        """
        octi_list = [(1000. * i, 1000. * i + 500.) for i in range(20)]
        pattern = xOnOrbitCalibrationPattern(octi_list, None)
        print(pattern)
        self.assertEqual(pattern.num_calibration_runs(1), 7)
        self.assertAlmostEqual(pattern.total_calibration_time(1), 3500.)
        self.assertEqual(pattern.num_calibration_runs(2), 7)
        self.assertAlmostEqual(pattern.total_calibration_time(2), 3500.)
        self.assertEqual(pattern.num_calibration_runs(3), 6)
        self.assertAlmostEqual(pattern.total_calibration_time(3), 3000.)



if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python
#
# Copyright (C) 2015--2018, the ixpeobssim team.
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
import datetime

from matplotlib.dates import date2num

from ixpeobssim.utils.time_ import MISSION_START_UNIX_TIME, current_met,\
    unix_to_string_utc, met_to_mjd, unix_to_met, met_to_unix, unix_to_string,\
    string_to_unix_utc, string_to_met_utc, current_time, string_to_unix_local,\
    unix_to_string_local, met_to_string_utc, string_to_met_local,\
    met_to_string_local, _format_datetime_lazy, met_to_num, MISSION_START_DATETIME,\
    xTimeInterval



class TestTime(unittest.TestCase):


    def test_intervals(self, start_met=0., stop_met=1000., padding=10.):
        """
        """
        dt1 = xTimeInterval(start_met, stop_met)
        print(dt1)

    def test_start(self):
        """Basic test.
        """
        t = MISSION_START_UNIX_TIME
        self.assertEqual(unix_to_met(t), 0)
        self.assertAlmostEqual(unix_to_met(float(t)), 0.)

    def test_start_datetime(self):
        """Test the start Unix time.
        """
        dt = (datetime.datetime(2017, 1, 1) - \
              datetime.datetime(1970, 1, 1)).total_seconds()
        self.assertEqual(dt, MISSION_START_UNIX_TIME)

    def test_met(self):
        """Test the current MET.
        """
        met1 = (datetime.datetime.utcnow() -\
                datetime.datetime(2017, 1, 1)).total_seconds()
        met2 = current_met()
        self.assertAlmostEqual(met1, met2, places=0)

    def test_datetime_fmt(self):
        """Test the conversion from time to a formatted date and time.
        """
        string = unix_to_string_utc(MISSION_START_UNIX_TIME)
        self.assertEqual(string, '2017-01-01T00:00:00.000000')

    def test_mjd(self):
        """Test the conversion to MJD.
        """
        self.assertAlmostEqual(met_to_mjd(0), 57754)

    def test_met_unix(self):
        """Convert back and forth from Unix time to MET.
        """
        t = current_time()
        met = unix_to_met(t)
        ut = met_to_unix(met)
        self.assertAlmostEqual(t, ut)

    def test_met_unix_string(self):
        """Convert back and forth from Unix time to MET strings.
        """
        ut = current_time()
        _ut = string_to_unix_utc(unix_to_string_utc(ut))
        self.assertAlmostEqual(_ut, ut, places=6)
        _ut = string_to_unix_local(unix_to_string_local(ut))
        self.assertAlmostEqual(_ut, ut, places=6)
        met = unix_to_met(ut)
        _met = string_to_met_utc(met_to_string_utc(met))
        self.assertAlmostEqual(_met, met, places=6)
        _met = string_to_met_local(met_to_string_local(met))
        self.assertAlmostEqual(_met, met, places=6)

    def test_string_to_time(self):
        """Test the conversion from string to timestamp.
        """
        string = unix_to_string(MISSION_START_UNIX_TIME)
        ut = string_to_unix_utc(string)
        met = string_to_met_utc(string)
        self.assertEqual(ut, MISSION_START_UNIX_TIME)
        self.assertEqual(met, 0)

    def test_lazy(self):
        """
        """
        target = '2017-01-01T00:00:00.0'
        self.assertEqual(_format_datetime_lazy(target), target)
        self.assertEqual(_format_datetime_lazy('2017-01-01'), target)
        self.assertEqual(_format_datetime_lazy('2017-01-01T00'), target)
        self.assertEqual(_format_datetime_lazy('2017-01-01T00:00'), target)
        self.assertEqual(_format_datetime_lazy('2017-01-01T00:00:00'), target)

    def test_met_to_num(self):
        """
        """
        self.assertEqual(met_to_num(0.), date2num(MISSION_START_DATETIME))



if __name__ == '__main__':
    unittest.main()

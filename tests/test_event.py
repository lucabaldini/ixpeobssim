#!/usr/bin/env python
#
# Copyright (C) 2019--2022, the ixpeobssim team.
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

import numpy
import unittest

from ixpeobssim.evt.event import xEventList
from ixpeobssim.utils.time_ import string_to_met


numpy.random.seed(666)


class TestEventTime(unittest.TestCase):

    LATE_MET = string_to_met('2025-01-01T00:00:00.0')

    """Unit test for the split of the event times in seconds and microseconds.
    """

    def _basic_test(self, t0, t1, size=1000000):
        """Basic test definition.
        """
        t = numpy.random.uniform(t0, t1, size)
        seconds, microseconds = xEventList.split_event_time(t)
        self.assertTrue((seconds >= t0).all())
        self.assertTrue((seconds <= t1).all())
        self.assertTrue((microseconds >= 0).all())
        self.assertTrue((microseconds <= 1000000).all())

    def test_basic(self):
        """
        """
        self._basic_test(0., 1.)

    def test_end_of_mission(self, duration=100000.):
        """
        """
        t0 = self.LATE_MET
        t1 = t0 + duration
        self._basic_test(t0, t1)

    def _test_round_trip(self, t0, t1, size=10):
        """
        """
        s = numpy.random.uniform(t0, t1, size).astype(numpy.int32)
        us = numpy.random.uniform(0, 999999.9, size).astype(numpy.int32)
        t = s + 1.e-6 * us
        # Add 100 ns to avoid numerical inaccuracies in the forward direction.
        t += 1.e-7
        seconds, microseconds = xEventList.split_event_time(t)
        self.assertTrue((s == seconds).all())
        self.assertTrue((us == microseconds).all())

    def test_precision(self, duration=100000.):
        """
        """
        t0 = self.LATE_MET
        t1 = t0 + duration
        self._test_round_trip(t0, t1)

    def test_sort_check(self):
        """Test the routine that check whether an array is sorted.
        """
        a = numpy.linspace(0., 100., 101)
        self.assertTrue(xEventList.array_is_sorted(a))
        b = numpy.random.uniform(size=10000)
        self.assertFalse(xEventList.array_is_sorted(b))

    def test_num_events(self):
        """Test all the bookkeeping handling the number of events.
        """
        size = 101
        time_ = numpy.linspace(0., 100., size)
        event_list = xEventList(time_)
        self.assertEqual(event_list.num_events(), size)



if __name__ == '__main__':
    unittest.main()

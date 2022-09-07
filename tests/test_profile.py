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


"""Unit test for the utils.profile module.
"""


import unittest
import time

from ixpeobssim.utils.profile import xChrono, xMemoryProfiler


class TestChrono(unittest.TestCase):

    """Unit test for the profile.xChrono class.
    """

    def test_basic(self, sleep_time=0.1):
        """Basic test of the interfaces.
        """
        c = xChrono()
        time.sleep(sleep_time)
        self.assertTrue(c() >= sleep_time)
        self.assertTrue(isinstance(str(c), str))

    def test_memory(self):
        """
        """
        prof = xMemoryProfiler()
        print(prof)



if __name__ == '__main__':
    unittest.main()

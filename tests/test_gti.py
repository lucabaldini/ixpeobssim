#!/usr/bin/env python
#
# Copyright (C) 2022, the ixpeobssim team.
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

from ixpeobssim.evt.gti import xGTIList


class TestGti(unittest.TestCase):

    def test_gti_complement(self):
        """
        """
        gti_list = xGTIList(0., 50., (0., 10.), (20., 30.), (40., 50.))
        print(gti_list)
        complement = gti_list.complement()
        print(complement)
        self.assertEqual(len(complement), len(gti_list) - 1)
        self.assertEqual(tuple(complement), ((10., 20.), (30., 40.)))



if __name__ == '__main__':
    unittest.main()

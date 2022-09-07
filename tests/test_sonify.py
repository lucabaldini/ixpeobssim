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

import numpy

from ixpeobssim.evt.sonify import xMidiNote, xMusicalScale
from ixpeobssim.utils.logging_ import logger


class TestSonify(unittest.TestCase):

    """Unit test for the sonification code.
    """

    def test_notes(self):
        """Test
        """
        for i in range(21, 128):
            note = xMidiNote(i)
            print(i, note)
        self.assertEqual(xMidiNote(21).name(), 'A0')
        self.assertEqual(xMidiNote(69).name(), 'A4')
        self.assertEqual(xMidiNote(127).name(), 'G9')

    def test_scales(self):
        """
        """
        for mode in xMusicalScale.SCALE_DICT:
            print(xMusicalScale(mode))




if __name__ == '__main__':
    unittest.main()

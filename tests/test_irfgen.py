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

from ixpeobssim.irf.ebounds import ENERGY_GRID, ENERGY_BINNING, NUM_CHANNELS, TLMIN, PI_BINNING
from ixpeobssim.irf.ebounds import channel_to_energy, energy_to_channel, digitize_channel


class TestIrfgen(unittest.TestCase):

    """Unit test for the irf __init__.py module.
    """

    def test_base(self):
        """
        """
        print(ENERGY_GRID)
        print(ENERGY_BINNING)
        print(PI_BINNING)
        self.assertEqual(NUM_CHANNELS, 375)
        self.assertEqual(TLMIN, 0)

    def test_conversion(self):
        """
        """
        energy = numpy.array([1.98000001, 2., 2.0199999])
        pi = energy_to_channel(energy)
        pi = digitize_channel(pi)
        print(energy)
        print(pi)
        print(channel_to_energy(pi))



if __name__ == '__main__':
    unittest.main()

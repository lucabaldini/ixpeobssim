#!/usr/bin/env python
#
# Copyright (C) 2018, the ixpeobssim team.
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


"""Unit test for the irf.mrf module.
"""

import unittest
import sys

from ixpeobssim.irf import load_arf, load_modf, load_mrf, DEFAULT_IRF_NAME
from ixpeobssim.irf.ebounds import ENERGY_GRID
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
if sys.flags.interactive:
    plt.ion()


class TestMrf(unittest.TestCase):

    """Unit test for the IXPE modulation response.
    """

    def base_test_mrf(self, irf_name, du_id=1):
        """
        """
        energy = ENERGY_GRID
        arf = load_arf(irf_name, du_id)
        modf = load_modf(irf_name, du_id)
        mrf = load_mrf(irf_name, du_id)
        delta = abs(arf(energy) * modf(energy) - mrf(energy)) / mrf(energy)
        plt.figure('Test MRF %s' % irf_name)
        plt.plot(energy, delta)
        self.assertTrue(delta.max() < 1e-5)

    def test_ixpe_default_mrf(self):
        """
        """
        self.base_test_mrf(DEFAULT_IRF_NAME)




if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

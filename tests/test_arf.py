#!/usr/bin/env python
#
# Copyright (C) 2016--2019, the ixpeobssim team.
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


"""Unit test for the irf.arf module.
"""

import unittest
import sys

import numpy

from ixpeobssim.instrument import DU_IDS
from ixpeobssim.irfgen.gpd import GPD_FILL_TEMPERATURE, quantum_efficiency
from ixpeobssim.irfgen.mma import MMA_AEFF_FILE_PATH
from ixpeobssim.irfgen.du import uv_filter_transparency
from ixpeobssim.irf import load_arf, DEFAULT_IRF_NAME
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()



class TestIxpeArf(unittest.TestCase):

    """Unit test for the IXPE effective area.
    """

    def setUp(cls):
        """Hack to prevent spurious warnings from numpy, see
        https://github.com/ContinuumIO/anaconda-issues/issues/6678

        also also, spurious import warning from Python itself.
        https://github.com/cython/cython/issues/1720
        """
        import warnings
        warnings.filterwarnings('ignore', message='numpy.dtype size changed')
        warnings.filterwarnings('ignore', message='can\'t resolve package')

    @unittest.skip('Combined IRF deprecated since v4')
    def test_one_to_three(self, tolerance=1.e-6):
        """Make sure the effective area for three DU is the sum of the effective
        areas for the three DUs.
        """
        aeff1 = load_arf(du_id=1).y
        aeff2 = load_arf(du_id=2).y
        aeff3 = load_arf(du_id=3).y
        aeff123 = load_arf(du_id=123).y
        delta = abs((aeff1 + aeff2 + aeff3) / aeff123 - 1.)
        self.assertTrue(numpy.allclose(delta, 0.0, rtol=tolerance, atol=tolerance))



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

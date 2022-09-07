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


"""Unit test for the energy dispersion facility.
"""

import numpy
import unittest
import sys

from ixpeobssim.instrument import DU_IDS
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.irfgen.gpd import _gpd_data_path
from ixpeobssim.irf import load_rmf, DEFAULT_IRF_NAME


if sys.flags.interactive:
    plt.ion()

class TestEnergyDispersion(unittest.TestCase):

    """Unit test for xModulationFactor.
    """

    @classmethod
    def setUpClass(cls):
        """Setup.

        Also, hack to prevent spurious warnings from numpy, see
        https://github.com/ContinuumIO/anaconda-issues/issues/6678

        also also, spurious import warning from Python itself.
        https://github.com/cython/cython/issues/1720
        """
        import warnings
        warnings.filterwarnings('ignore', message='numpy.dtype size changed')
        warnings.filterwarnings('ignore', message='can\'t resolve package')
        cls.interactive = sys.flags.interactive

    def test_rvs(self, num_events=100000):
        """
        """
        mc_energy = 5.
        _e = numpy.full(num_events, mc_energy)
        _bins = numpy.linspace(0, 256, 257) - 0.5
        for du_id in DU_IDS:
            edisp = load_rmf(DEFAULT_IRF_NAME, du_id)
            _slice = edisp.matrix.slice(mc_energy)
            _ch = edisp.matrix.rvs(_e)
            plt.figure()
            n, bins, patches = plt.hist(_ch, bins=_bins, histtype='step')
            bin_width = (bins[1] - bins[0])
            scale = n.sum()*bin_width/_slice.norm()
            plt.plot(_slice.x, scale*_slice.y)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

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

"""Unit test for the irf.rmf module.
"""

from __future__ import print_function, division

import unittest
import sys

import numpy

from ixpeobssim.instrument import DU_IDS
from ixpeobssim.irfgen.gpd import _gpd_data_path
from ixpeobssim.irf import load_rmf, DEFAULT_IRF_NAME
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()



class TestIxpeRmf(unittest.TestCase):

    """Unit test for the IXPE energy dispersion.
    """

    @classmethod
    def setUpClass(cls):
        """Setup. Create a few objects to be used for testing.

        Also, hack to prevent spurious import warning from Python itself.
        https://github.com/cython/cython/issues/1720
        """
        import warnings
        warnings.filterwarnings('ignore', message='can\'t resolve package')

    def test_plot(self, du_id=1):
        """Plot the energy dispersion.
        """
        edisp = load_rmf(DEFAULT_IRF_NAME, du_id)
        edisp.plot()

    def test_ixpe_rmf_matrix_norm(self):
        """Test the IXPE energy dispersion normalization.

        Take a number of vertical slices of the energy dispersion matrix
        and make sure they are normalized to unity.
        """
        for du_id in DU_IDS:
            edisp = load_rmf(DEFAULT_IRF_NAME, du_id)
            emin = edisp.matrix.ymin()
            emax = edisp.matrix.ymax()
            de = emax - emin
            for energy in numpy.linspace(emin + 0.25 * de, emax - 0.25 * de, 10):
                _delta = abs(edisp.matrix.hslice(energy).norm() - 1)
                self.assertTrue(_delta < 1e-3, 'diff. %.9f in DU%d' %\
                                (_delta, du_id))

    def test_conversions(self):
        """Test some basic channel <-> energy conversions.
        """
        edisp = load_rmf(DEFAULT_IRF_NAME, du_id=1)
        chan = numpy.arange(100)
        energy = edisp.channel_to_energy(chan)
        _chan = edisp.energy_to_channel(energy)
        self.assertTrue(numpy.allclose(chan, _chan))
        delta = 0.499999 * (energy[1] - energy[0])
        _chan = edisp.energy_to_channel(energy + delta)
        self.assertTrue(numpy.allclose(chan, _chan))
        _chan = edisp.energy_to_channel(energy - delta)
        self.assertTrue(numpy.allclose(chan, _chan))



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

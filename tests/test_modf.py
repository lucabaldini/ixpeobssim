#!/usr/bin/env python
#
# Copyright (C) 2016--2018, the ixpeobssim team.
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


"""Unit test for the irf.modf module.
"""

import unittest
import sys

import numpy

from ixpeobssim.instrument import DU_IDS
from ixpeobssim.irfgen.gpd import _gpd_data_path
from ixpeobssim.irf import load_modf, DEFAULT_IRF_NAME
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()

GPD_MODF_FILE_PATH = _gpd_data_path('modf_hedme8020_1atm_1cm.txt')


class TestIxpeModf(unittest.TestCase):

    """Unit test for the IXPE modulation factor.
    """

    @classmethod
    def setUpClass(cls):
        """Setup.
        """
        cls.modf_dict = {}
        for du_id in DU_IDS:
            cls.modf_dict[du_id] = load_modf(DEFAULT_IRF_NAME, du_id)

    def test_plot(self):
        """Basic plotting test.
        """
        for du_id in DU_IDS:
            self.modf_dict[du_id].plot()

    @unittest.skip('only applicable to IRFs v2 and before')
    def test_modf(self):
        """Test the IXPE modulation factor.

        This loads the modulation factor from the IXPE modf FITS file,
        then the corresponding values from the text file that the response
        functions are created from and makes sure that the two are close
        enough.
        """
        _x, _y = numpy.loadtxt(GPD_MODF_FILE_PATH, unpack=True)
        _mask = _y > 0.
        _x = _x[_mask]
        _y = _y[_mask]
        for du_id in DU_IDS:
            _delta = abs((_y - self.modf_dict[du_id](_x))/_y)
            self.assertTrue(_delta.max() < 7e-3, 'max. diff. %.9f in DU%d' %\
                            (_delta.max(), du_id))

    def test_unphysical_polarization(self, du_id=1):
        """Make sure that we handle correctly cases where we pass a polarization
        degree smaller than zero or greater than one.
        """
        modf = self.modf_dict[du_id]
        energy = numpy.linspace(2., 8., 7)
        pd = (energy - 5.) / 2.
        pa = 0.
        self.assertRaises(SystemExit, modf.rvs_phi, energy, pd, pa)
        pd = energy / 5.
        self.assertRaises(SystemExit, modf.rvs_phi, energy, pd, pa)

    def test_extrapolation(self):
        """Test that the modulation factor is positive in the entire 0--15 keV
        energy band.
        """
        energy = numpy.linspace(0., 15., 1000)
        for du_id in DU_IDS:
            modf = self.modf_dict[du_id](energy)
            self.assertTrue((modf >= 0.).all())


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

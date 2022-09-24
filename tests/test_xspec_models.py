#!/usr/bin/env python
#
# Copyright (C) 2015--2019, the ixpeobssim team.
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


"""Unit test for the srcmodel.spectrum.xXspecModel class.
"""

import unittest
import os
import sys

import numpy
import matplotlib.pyplot as plt

from ixpeobssim import IXPEOBSSIM_SRCMODEL, IXPEOBSSIM_DATA
from ixpeobssim.utils.environment import PYXSPEC_INSTALLED
from ixpeobssim.srcmodel.spectrum import power_law, cutoff_power_law, xXspecModel
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
if sys.flags.interactive:
    plt.ion()


class TestXspecModels(unittest.TestCase):

    """Unit test for xXspecModel.
    """

    @classmethod
    def setUpClass(cls):
        """Setup.
           Create a few objects to be used for testing.
        """
        if not PYXSPEC_INSTALLED:
            return
        cls.pl_index = 1.
        cls.pl_norm = 10.
        cls.spec = xXspecModel('powerlaw', [cls.pl_index, cls.pl_norm])

    @unittest.skip('Fix issue#382')
    def test_power_law(self):
        """Test if xXspecModel reproduces the right powerlaw distribution.
        """
        if not PYXSPEC_INSTALLED:
            return
        pl = power_law(self.pl_norm, self.pl_index)
        args = 'powerlaw', [self.pl_index, self.pl_norm]
        spec_corr = xXspecModel(*args, correct=True)
        spec_raw = xXspecModel(*args, correct=False)
        plt.figure('XSPEC powerlaw accuracy')
        x = spec_corr.x
        y = spec_corr.y
        _delta = abs(pl(x) - y) / y
        self.assertTrue(_delta.max() < 1.e-5, 'max. diff. %.4f' % _delta.max())
        plt.plot(x, _delta, label='Corrected')
        y = spec_raw.y
        _delta = abs(pl(x) - y) / y
        plt.plot(x, _delta, label='Uncorrected')
        setup_gca(xlabel='Energy [keV]', ylabel='Fractional deviation',
                  logx=True, logy=True, grids=True, xmin=1., xmax=12., legend=True)

    @unittest.skip('Fix issue#382')
    def test_cutoff_power_law(self, cutoff_energy=4.):
        """Test if xXspecModel reproduces the right powerlaw distribution.
        """
        if not PYXSPEC_INSTALLED:
            return
        pl = cutoff_power_law(self.pl_norm, self.pl_index, cutoff_energy)
        args = 'cutoffpl', [self.pl_index, cutoff_energy, self.pl_norm]
        spec_corr = xXspecModel(*args, correct=True)
        spec_raw = xXspecModel(*args, correct=False)
        plt.figure('XSPEC cutoff powerlaw accuracy')
        x = spec_corr.x
        y = spec_corr.y
        _delta = abs(pl(x) - y) / y
        self.assertTrue(_delta.max() < 1.e-5, 'max. diff. %.4f' % _delta.max())
        plt.plot(x, _delta, label='Corrected')
        y = spec_raw.y
        _delta = abs(pl(x) - y) / y
        plt.plot(x, _delta, label='Uncorrected')
        setup_gca(xlabel='Energy [keV]', ylabel='Fractional deviation',
                  logx=True, logy=True, grids=True, xmin=1., xmax=12., legend=True)

    def _test_param_dict(self):
        """Test initializing the parameters from a dictionary rather than a list.
        """
        if not PYXSPEC_INSTALLED:
            return
        spec = xXspecModel('powerlaw', {1: self.pl_index, 2: self.pl_norm})
        self.assertTrue(numpy.allclose(spec.y, self.spec.y))

    def test_read_parfile(self):
        """Test xspec model buiding method.
        """
        if not PYXSPEC_INSTALLED:
            return
        file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'ascii', 'powerlaw_parfile.par')
        model = xXspecModel.from_file(file_path)
        print(model)

    def test_save(self):
        """Test saving a model.
        """
        if not PYXSPEC_INSTALLED:
            return
        file_path = os.path.join(IXPEOBSSIM_DATA, 'test_xspec_model.txt')
        self.spec.save_ascii(file_path)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

#!/usr/bin/env python
#
# Copyright (C) 2023, the ixpeobssim team.
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


"""Unit test for the irf.arf module in the grayfilter version.
"""

import unittest
import sys

import numpy

from ixpeobssim.instrument import DU_IDS
from ixpeobssim.irf import load_arf
from ixpeobssim.utils.matplotlib_ import plt, setup_gca

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

    def test_aeff(self):
        """Test the new padding for the spectral response files introduced
        after https://github.com/lucabaldini/ixpeobssim/issues/713
        """
        for du_id in DU_IDS:
            aeff = load_arf('ixpe:obssim:v12', du_id, gray_filter=True)
            plt.figure('Effective area with the gray filter (DU %d)' % du_id)
            plt.plot(aeff.x, aeff.y, 'o')
            E = numpy.linspace(1., 12., 1000)
            plt.plot(E, aeff(E))
            setup_gca(xlabel='Energy [keV]', ylabel='Effective area [cm$^2$]', logx=True,
                logy=True)


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

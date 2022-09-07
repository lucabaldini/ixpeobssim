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


import unittest
import sys

import numpy

from ixpeobssim.irf import load_arf, load_rmf, DEFAULT_IRF_NAME
from ixpeobssim.irf.ebounds import energy_to_channel
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color

if sys.flags.interactive:
    plt.ion()



class TestEnergyRedistribution(unittest.TestCase):

    """Unit test for the IXPE energy redistribution.
    """

    def test(self, irf_name=DEFAULT_IRF_NAME):
        """
        """
        aeff = load_arf(irf_name)
        edisp = load_rmf(irf_name)
        plt.figure('Energy dispersion matrix')
        edisp.matrix.plot()
        plt.figure('Energy dispersion direct')
        for energy in (2., 4., 6., 9.):
            s = edisp.matrix.hslice(energy)
            s.plot(label='%.1f keV' % energy)
        setup_gca(legend=True, grids=True, logy=True, ymin=1.e-4, ylabel='Probability density')
        plt.figure('Energy dispersion inverse')
        for energy in (2., 4., 6., 9.):
            channel = energy_to_channel(energy)
            s = edisp.matrix.vslice(channel)
            x = s.x
            y = s(x)
            y /= y.sum()
            plt.plot(x, y, label='%.1f keV' % energy)
            y = s(x) * aeff(x) * x**-2.
            y /= y.sum()
            plt.plot(x, y, color=last_line_color(), ls='dashed')
        setup_gca(legend=True, grids=True, logy=True, ymin=1.e-6,
                  xlabel='Energy [keV]', ylabel='Probability density')


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

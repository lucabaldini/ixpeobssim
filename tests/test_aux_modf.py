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
import sys

from ixpeobssim.irfgen.gpd import xModfDataInterface
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color, residual_plot

if sys.flags.interactive:
    plt.ion()


class TestAuxModf(unittest.TestCase):

    """
    """

    @classmethod
    def setUpClass(cls):
        """Load the modulation factor data.
        """
        cls.modf_data = xModfDataInterface()
        plt.figure('Modulation factor spline')
        cls.modf_data.plot()

    def test_raw_data(self):
        """
        """
        for i, pressure in enumerate(self.modf_data.pressure):
            ax1, ax2 = residual_plot('Modulation factor @ %d mbar' % pressure)
            s = self.modf_data.pressure_slice(pressure)
            s.plot(label='%d mbar' % pressure)
            x = self.modf_data.energy
            y = self.modf_data.mu[:,i]
            dy = self.modf_data.mu_err[:,i]
            plt.errorbar(x, y, dy, fmt='o', color=last_line_color())
            setup_gca(ylabel='Modulation factor', grids=True)
            plt.sca(ax2)
            plt.errorbar(x, y - s(x), dy, fmt='o')
            setup_gca(xlabel='Energy [keV]', ymin=-0.015, ymax=0.015, grids=True)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

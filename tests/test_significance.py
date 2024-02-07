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

import unittest
import sys

import numpy
import scipy.stats

from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import residual_plot, plt, setup_gca

if sys.flags.interactive:
    plt.ion()



class TestSignificance(unittest.TestCase):

    """Unit test for the significance calculation in kislat2015.

    See https://github.com/lucabaldini/ixpeobssim/issues/709 for more
    background information.
    """

    @staticmethod
    def _significance_old(confidence):
        """Old routine.
        """
        return scipy.stats.norm.ppf(confidence)

    @staticmethod
    def _significance_new(confidence):
        """New routine proposed by Phil.
        """
        return scipy.stats.norm.ppf(0.5 * confidence + 0.5)

    def test(self):
        """Quick test comparing the two.
        """
        ax1, ax2 = residual_plot('Confidence calculation')
        conf = numpy.linspace(0., 1., 1000)
        sig_old = self._significance_old(conf)
        sig_new = self._significance_new(conf)
        plt.plot(conf, sig_old, label='Old')
        plt.plot(conf, sig_new, label='New')
        setup_gca(xlabel='Confidence', ylabel='Significance', grids=True, legend=True)
        plt.sca(ax2)
        mask = conf > 0.6
        plt.plot(conf[mask], sig_new[mask] / sig_old[mask])
        plt.axhline(1., ls='dashed')
        setup_gca(xlabel='Confidence', ylabel='Ratio', xmin=0., ymin=0, ymax=2.5, grids=True)
        for i, p in enumerate((0, 0.682689492137086, 0.954499736103642, 0.997300203936740)):
            self.assertAlmostEqual(self._significance_new(p), i)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

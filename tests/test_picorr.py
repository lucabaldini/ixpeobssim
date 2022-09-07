#!/usr/bin/env python
#
# Copyright (C) 2022, the ixpeobssim team.
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

"""Unit tests for the PI scaling
"""


from __future__ import print_function, division

import sys
import unittest

import numpy

from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.evt.picorr import xTrapezoidPdf
from ixpeobssim.irf.ebounds import energy_to_channel
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, residual_plot

if sys.flags.interactive:
    plt.ion()

# pylint: disable=invalid-name


class TestPulseInvariantCorrection(unittest.TestCase):

    """Unit test for the pulse invariant correction.
    """

    @staticmethod
    def correct_pi(pi, scale, offset=0.):
        """
        """
        pi = pi + numpy.random.uniform(-0.5, 0.5, size=pi.shape)
        pi = pi * scale + offset
        pi = numpy.rint(pi)
        return pi

    def _test_scaling_base(self, energy, label, scale=1.2):
        """Base test function.
        """
        pi = energy_to_channel(energy)
        pi_scaled1 = energy_to_channel(energy * scale)
        pi_scaled2 = self.correct_pi(pi, scale)
        binning = numpy.linspace(-0.5, 374.5, 376)
        h1 = xHistogram1d(binning).fill(pi_scaled1)
        h2 = xHistogram1d(binning).fill(pi_scaled2)
        ax1, ax2 = residual_plot('%s scaling' % label)
        h1.plot(label='Energy scaling')
        h2.plot(label='PI scaling')
        setup_gca(xlabel='Channel', legend=True)
        mask = h1.content > 0
        delta = (h1.content[mask] - h2.content[mask]) / numpy.sqrt(h1.content[mask])
        chisq = (delta**2.).sum()
        dof = len(delta)
        logger.info('Chisquare: %.3f / %d dof', chisq, dof)
        plt.sca(ax2)
        plt.plot(h2.content / h1.content)
        setup_gca(xmin=-0.5, xmax=374.5, grids=True)

    def test_constant(self, emin=2., emax=8., num_events=1000000):
        """
        """
        energy = numpy.random.uniform(emin, emax, size=num_events)
        self._test_scaling_base(energy, 'constant')

    def test_triangular(self, emin=2., emax=3., num_events=1000000):
        """
        """
        emean = 0.5 * (emin + emax)
        energy = numpy.random.triangular(emin, emean, emax, size=num_events)
        self._test_scaling_base(energy, 'triangular')

    def test_gaussian(self, mean=3., sigma=0.25, num_events=1000000):
        """
        """
        energy = numpy.random.normal(mean, sigma, size=num_events)
        self._test_scaling_base(energy, 'gaussian')




if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

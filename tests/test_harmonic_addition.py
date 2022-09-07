#!/usr/bin/env python
#
# Copyright (C) 2015, the ixpeobssim team.
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
import numpy
import sys

from ixpeobssim.srcmodel.polarization import harmonic_addition
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.utils.logging_ import logger

if sys.flags.interactive:
    plt.ion()


def _pdf(phi, F=1., m=0.5, phi0=0):
    """Basic modulation curve definition (with arbitrary normalization).
    """
    return F * (1 + m * numpy.cos(2. * (phi - phi0))) / (2. * numpy.pi)


def _pdf_sum(phi, *params):
    """Brute-force function addition for an arbitrary number of components.

    (Using the same signature as ixpeobssim.polarization.harmonic_addition)
    """
    return sum(_pdf(phi, *p) for p in params)
    



class TestHarmonicAddition(unittest.TestCase):

    """Test class for the harmonic addition theorem.
    """

    def basic_twofold_test(self, F1, m1, delta1, F2, m2, delta2):
        """Basic harmonic addition test with two components.
        """
        logger.info('Testing harmonic addition...')
        logger.info('(%.4f, %.4f, %.4f) +' % (F1, m1, delta1))
        logger.info('(%.4f, %.4f, %.4f) =' % (F2, m2, delta2))
        phi = numpy.linspace(-numpy.pi, numpy.pi, 100)
        plt.figure()
        y1 = _pdf_sum(phi, (F1, m1, delta1), (F2, m2, delta2))
        F, m, delta = harmonic_addition((F1, m1, delta1), (F2, m2, delta2))
        logger.info('------------------------')
        logger.info('(%.4f, %.4f, %.4f)' % (F, m, delta))
        y2 = _pdf(phi, F, m, delta)
        self.assertTrue(max(abs(y2 - y1)) < 1e-12)
        plt.plot(phi, y1, label='Brute force')
        plt.plot(phi, y2, label='Harmonic addition')
        plt.axis([-numpy.pi, numpy.pi, 0, None])
        plt.legend()

    def test_points(self):
        """Basic test.
        """
        self.basic_twofold_test(0.5, 0.50, 0., 0.5, 0.50, 0.50 * numpy.pi)
        self.basic_twofold_test(0.5, 0.50, 0., 0.5, 0.25, 0.25 * numpy.pi)
        self.basic_twofold_test(0.5, 0.25, 0., 0.5, 0.50, 0.25 * numpy.pi)
        self.basic_twofold_test(3.5, 0.50, 0., 5.0, 0.20, 0.25 * numpy.pi)
        self.basic_twofold_test(1.5, 0.50, 0., 8.0, 0.20, 0.50 * numpy.pi)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

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

"""Unit tests for srcmodel.spectrum
"""

from __future__ import print_function, division

import sys
import unittest

import numpy

from ixpeobssim.srcmodel.spectrum import powerlaw, cutoffpl
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()


class TestModels(unittest.TestCase):

    """Unit tests for the spectral models.
    """

    @staticmethod
    def basic_test(model, figure_name):
        """
        """
        E = numpy.linspace(1., 10., 250)
        t = numpy.linspace(0., 1000., 3)
        fmt = dict(xlabel='Energy [keV]',
                   ylabel='dN/dE [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]')
        plt.figure(figure_name)
        for _t in t:
            s = xInterpolatedUnivariateSpline(E, model(E, _t), **fmt)
            s.plot(logx=True, logy=True, label='t = %.1f' % _t)
        plt.legend()

    @staticmethod
    def norm(t):
        """Time-dependent normalization.
        """
        return 1. + t * 0.001

    @staticmethod
    def index(t):
        """Time-dependent index.
        """
        return 2. - t * 0.001

    def test_powerlaw_stationary(self):
        """
        """
        model = powerlaw(self.norm(0), self.index(0))
        self.basic_test(model, 'Power law (stationary)')
        # Make sure in the stationary case we can call the model without passing
        # the time explicitly,
        E = numpy.linspace(1., 10., 10)
        self.assertTrue((model(E) == model(E, 0.)).all())


    def test_powerlaw_norm(self):
        """
        """
        model = powerlaw(self.norm, self.index(0))
        self.basic_test(model, 'Power law (time-dependent normalization)')

    def test_powerlaw_index(self):
        """
        """
        model = powerlaw(self.norm(0), self.index)
        self.basic_test(model, 'Power law (time-dependent index)')

    def test_powerlaw(self):
        """
        """
        model = powerlaw(self.norm, self.index)
        self.basic_test(model, 'Power law (time-dependent)')

    def test_cutoffpl(self):
        model = cutoffpl(self.norm, self.index, 5.)
        self.basic_test(model, 'Power law with exponential cutoff')





if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

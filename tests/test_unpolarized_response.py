#!/urs/bin/env python
#
# Copyright (C) 2019, the ixpeobssim team.
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

from ixpeobssim.irf.modf import xAzimuthalResponseGenerator
from ixpeobssim.core.modeling import xModulationCurveRad
from ixpeobssim.core.fitting import fit_histogram
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.utils.matplotlib_ import plt

import sys
if sys.flags.interactive:
    plt.ion()


"""We explictely set the random seed to have reproducible results.
"""
numpy.random.seed(0)


class TestUnpolarizedResponse(unittest.TestCase):

    """This unit test was triggered by this issue
    https://bitbucket.org/ixpesw/ixpeobssim/issues/182
    and is meant to make sure that, when generating unpolarized photons,
    the modulation curve looks statistically flat.
    """

    @classmethod
    def setUpClass(cls):
        """Setup.
        """
        cls.generator = xAzimuthalResponseGenerator()

    def test_generator(self, size=20000000, phase=0.):
        """
        """
        modulation = numpy.full(size, 0.)
        phi = self.generator.rvs_phi(modulation, phase)
        binning = numpy.linspace(-numpy.pi, numpy.pi, 100)
        hist = xHistogram1d(binning).fill(phi)
        model = xModulationCurveRad()
        fit_histogram(model, hist)
        if sys.flags.interactive:
            plt.figure()
            hist.plot()
            model.plot()
            model.stat_box()
        m = model.parameter_value('Modulation')
        dm = model.parameter_error('Modulation')
        delta = abs(m / dm)
        self.assertTrue(delta < 3.)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

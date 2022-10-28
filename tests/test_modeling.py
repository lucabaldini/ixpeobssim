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


import itertools
import unittest
import sys

from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.core.modeling import xConstant, xLine, xGaussian, xFe55,\
    xPowerLaw, xPowerLawExpCutoff, xPixelPha, xModulationCurveRad,\
    xModulationCurveDeg, xExponential, xExponentialOffset, xFitModelBase

if sys.flags.interactive:
    plt.ion()


class TestModeling(unittest.TestCase):

    """
    """

    def test_models(self):
        """Test all the models in the modeling module.
        """
        for model in [xConstant(), xLine(), xGaussian(), xFe55(),
                      xPowerLaw(), xPowerLawExpCutoff(), xPixelPha(),
                      xModulationCurveRad(), xModulationCurveDeg(),
                      xExponential(), xExponentialOffset()]:
            plt.figure()
            model.plot()
            model.stat_box()

    def test_model_addition(self):
        """Test addition of the modules
        """
        models = [xConstant(), xLine(), xGaussian(), xFe55(),
                  xPowerLaw(), xPowerLawExpCutoff(), xPixelPha(),
                  xModulationCurveRad(), xModulationCurveDeg(),
                  xExponential(), xExponentialOffset()]
        for first, second in itertools.combinations(models, 2):
            model = first + second
            self.assertEqual(len(model.bounds[0]), len(model.parameters))
            self.assertEqual(len(model.bounds[1]), len(model.parameters))
            if first.bounds != xFitModelBase.PARAMETER_DEFAULT_BOUNDS:
                self.assertEqual(model.bounds[0][:len(first.parameters)],
                                 first.bounds[0])
                self.assertEqual(model.bounds[1][:len(first.parameters)],
                                 first.bounds[1])
            if second.bounds != xFitModelBase.PARAMETER_DEFAULT_BOUNDS:
                self.assertEqual(model.bounds[0][len(first.parameters):],
                                 second.bounds[0])
                self.assertEqual(model.bounds[1][len(first.parameters):],
                                 second.bounds[1])

if __name__ == "__main__":
    unittest.main(exit=False)

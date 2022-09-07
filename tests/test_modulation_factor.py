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


"""Unit test for the modulation factor facility.
"""

import numpy
import unittest
import sys

from ixpeobssim import IXPEOBSSIM_DOC_FIGURES
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.irfgen.gpd import _gpd_data_path
from ixpeobssim.irf import load_modf, DEFAULT_IRF_NAME
from ixpeobssim.core.modeling import xModulationCurveRad
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.utils.misc import pairwise_enum

if sys.flags.interactive:
    plt.ion()


GPD_MODF_FILE_PATH = _gpd_data_path('modf_hedme8020_1atm_1cm.txt')


"""We explictely set the random seed to have reproducible results.
"""
numpy.random.seed(0)


class TestModulationFactor(unittest.TestCase):

    """Unit test for xModulationFactor.
    """

    @classmethod
    def setUpClass(cls, du_id=1):
        """Setup---here we essentially create the modulation factor.

        Also, hack to prevent spurious import warning from Python itself.
        https://github.com/cython/cython/issues/1720
        """
        import warnings
        warnings.filterwarnings('ignore', message='can\'t resolve package')
        cls.modf = load_modf(DEFAULT_IRF_NAME, du_id)

    def test_constant(self, num_events=1000000, polarization_degree=1.,
                      polarization_angle=numpy.radians(20.)):
        """Test the modulation factor as a random number generator when
        both the polarization angle and degrees are energy- and
        time-independent.
        """
        poldegree = numpy.full(num_events, polarization_degree)
        polangle = numpy.full(num_events, polarization_angle)
        self.modf.generator.plot()

        emin = self.modf.xmin()
        emax = self.modf.xmax()
        energy = numpy.random.uniform(emin, emax, num_events)
        phi = self.modf.rvs_phi(energy, poldegree, polangle)
        ebinning = numpy.linspace(emin, emax, 10)
        phi_binning = numpy.linspace(-numpy.pi, numpy.pi, 100)
        fit_models = []
        for i, (_emin, _emax) in pairwise_enum(ebinning):
            plt.figure()
            _emean = 0.5*(_emin + _emax)
            _mask = (energy > _emin)*(energy < _emax)
            _phi = phi[_mask]
            _hist = xHistogram1d(phi_binning).fill(_phi)
            _model = xModulationCurveRad()
            _hist.fit(_model)
            _model.emean = _emean
            fit_models.append(_model)
            _hist.plot()
            _model.plot(label='Energy: %.2f--%.2f keV' % (_emin, _emax))

        plt.figure()
        _x = [_model.emean for _model in fit_models]
        _y = [_model.parameter_value('Phase') for _model in fit_models]
        _dy = [_model.parameter_error('Phase') for _model in fit_models]
        plt.errorbar(_x, _y, yerr=_dy, fmt='o')
        plt.plot(_x, numpy.array([polarization_angle]*len(_x)))
        plt.xlabel('Energy [keV]')
        plt.ylabel('Modulation angle [$^\circ$]')

        plt.figure()
        _y = [_model.parameter_value('Modulation') for _model in fit_models]
        _dy = [_model.parameter_error('Modulation') for _model in fit_models]
        plt.errorbar(_x, _y, yerr=_dy, fmt='o')
        plt.axis([emin, emax, 0, 1])
        self.modf.plot()
        plt.xlabel('Energy [keV]')
        plt.ylabel('Modulation modulation')



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

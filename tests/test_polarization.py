#!/urs/bin/env python
#
# Copyright (C) 2015--2019, the ixpeobssim team.
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

from ixpeobssim.irf.modf import xAzimuthalResponseGenerator
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.irf import load_arf, load_modf, DEFAULT_IRF_NAME
from ixpeobssim.core.rand import xUnivariateGenerator
from ixpeobssim.core.modeling import xModulationCurveRad
from ixpeobssim.core.fitting import fit_histogram
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.srcmodel.polarization import fourier_series_factory,\
    xRadialPolarizationField, xTangentialPolarizationField
from ixpeobssim.utils.logging_ import logger

import sys
if sys.flags.interactive:
    plt.ion()


"""
"""

class TestPolarization(unittest.TestCase):

    """
    """

    @classmethod
    def setUpClass(cls, du_id=1):
        """Setup.
        """
        cls.aeff = load_arf(DEFAULT_IRF_NAME, du_id)
        cls.modf = load_modf(DEFAULT_IRF_NAME, du_id)
        cls.generator = xAzimuthalResponseGenerator()
        cls.emin = 1.
        cls.emax = 10.

    def polarization_degree(self, energy):
        """
        """
        return (1. - (energy - self.emin)/(self.emax - self.emin))

    def test_simplest(self, size=100000, phase=0.):
        """
        """
        modf = 1.
        energy = numpy.random.uniform(self.emin, self.emax, size)
        modulation = modf*self.polarization_degree(energy)
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
        fit_degree = model.parameter_value('Modulation')
        fit_error = model.parameter_error('Modulation')
        mean_energy = numpy.mean(energy)
        exp_degree = self.polarization_degree(mean_energy)
        delta = abs(fit_degree - exp_degree)/fit_error
        msg = 'delta = %.3f' % delta
        self.assertTrue(delta < 5., msg)

    def test_uniform(self, size=100000, phase=0.):
        """
        """
        energy = numpy.random.uniform(self.emin, self.emax, size)
        modulation = self.modf(energy)*self.polarization_degree(energy)
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
        mu_effective = self.modf.weighted_average(energy)
        fit_degree = model.parameter_value('Modulation') / mu_effective
        fit_error = model.parameter_error('Modulation') / mu_effective
        exp_degree = (self.polarization_degree(energy)*self.modf(energy))\
                     .sum()/self.modf(energy).sum()
        delta = abs(fit_degree - exp_degree)/fit_error
        msg = 'delta = %.3f' % delta
        self.assertTrue(delta < 5., msg)

    def test_power_law(self, size=1000000, phase=0., index=2.):
        """
        """
        _x = numpy.linspace(self.emin, self.emax, 100)
        _y = _x**(-index)
        generator = xUnivariateGenerator(_x, _y)
        energy = generator.rvs(size)
        modulation = self.modf(energy)*self.polarization_degree(energy)
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
        mu_effective = self.modf.weighted_average(energy)
        fit_degree = model.parameter_value('Modulation') / mu_effective
        fit_error = model.parameter_error('Modulation') / mu_effective
        exp_degree = (self.polarization_degree(energy)*self.modf(energy))\
                     .sum()/self.modf(energy).sum()
        delta = abs(fit_degree - exp_degree)/fit_error
        msg = 'delta = %.3f' % delta
        self.assertTrue(delta < 5., msg)

    def test_exponential_law(self, size=1000000, phase=0., base=10.):
        """
        """
        _x = numpy.linspace(self.emin, self.emax, 100)
        _y = base**(-_x)
        generator = xUnivariateGenerator(_x, _y)
        energy = generator.rvs(size)
        modulation = self.modf(energy)*self.polarization_degree(energy)
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
        mu_effective = self.modf.weighted_average(energy)
        fit_degree = model.parameter_value('Modulation') / mu_effective
        fit_error = model.parameter_error('Modulation') / mu_effective
        exp_degree = (self.polarization_degree(energy)*self.modf(energy))\
                     .sum()/self.modf(energy).sum()
        delta = abs(fit_degree - exp_degree)/fit_error
        msg = 'delta = %.3f' % delta
        self.assertTrue(delta < 5., msg)

    def test_fourier_expansion(self):
        """
        """
        factory = fourier_series_factory((0.2, 0.5), (0.3, 0.5))
        x = numpy.linspace(0., 1., 100)
        y1 = factory(x)
        y2 = 1. + 0.2 * numpy.cos(2 * numpy.pi * (x - 0.5)) +\
             0.3 * numpy.cos(2 * numpy.pi * 2. * (x - 0.5))
        self.assertFalse((y2 - y1).any())
        if sys.flags.interactive:
            plt.figure()
            plt.plot(x, y1)

    def test_radial_field(self):
        """Test the basic xRadialPolarizationField functionality.
        """
        ra0, dec0 = 30., 30.
        delta = 0.1
        epsilon = 1.e-12
        field = xRadialPolarizationField(ra0, dec0)
        for args, target in (
            ((ra0, dec0 + delta), 0.),
            ((ra0 + delta, dec0), 0.5 * numpy.pi),
            ((ra0, dec0 - delta), numpy.pi),
            ((ra0 - delta, dec0), -0.5 * numpy.pi),
            ((ra0 - epsilon, dec0 - delta), -numpy.pi)
            ):
            logger.info('Testing (%s, %s) -> %s', *args, target)
            self.assertAlmostEqual(field.polarization_angle(*args), target)

    def test_tangential_field(self):
        """Test the basic xTangentialPolarizationField functionality.
        """
        ra0, dec0 = 30., 30.
        delta = 0.1
        epsilon = 1.e-12
        field = xTangentialPolarizationField(ra0, dec0)
        for args, target in (
            ((ra0, dec0 + delta), 0.5 * numpy.pi),
            ((ra0 + delta, dec0), numpy.pi),
            ((ra0, dec0 - delta), -0.5 * numpy.pi),
            ((ra0 - delta, dec0), 0.),
            ((ra0 + delta, dec0 - epsilon), -numpy.pi)
            ):
            self.assertAlmostEqual(field.polarization_angle(*args), target)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

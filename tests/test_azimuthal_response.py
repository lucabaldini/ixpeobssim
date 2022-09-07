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


"""Unit test for the count spectrum facility.
"""

import numpy
import unittest
import sys

from ixpeobssim import IXPEOBSSIM_DOC_FIGURES
from ixpeobssim.irf.modf import xAzimuthalResponseGenerator
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.core.modeling import xModulationCurveRad
from ixpeobssim.core.fitting import fit_histogram
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.utils.logging_ import logger

if sys.flags.interactive:
    plt.ion()


"""We explictely set the random seed to have reproducible results.
"""
numpy.random.seed(0)


class TestAzimuthalResponse(unittest.TestCase):

    """Unit test for xAzimuthalResponseGenerator
    """

    @classmethod
    def setUpClass(cls):
        """Setup---here we essentially create the effective area.
        """
        cls.generator = xAzimuthalResponseGenerator()

    def test_basics(self):
        """Plot the basic stuff.
        """
        plt.figure('Azimuthal generator pdf')
        self.generator.plot()
        plt.figure('Azimuthal generator ppf')
        self.generator.ppf.plot()

    def test_norm(self):
        """Test the normalization.
        """
        for m in [0, 0.5, 1]:
            _x = numpy.linspace(0, 2 * numpy.pi, 100)
            _y = self.generator.pdf(_x, m)
            s = xInterpolatedUnivariateSplineLinear(_x, _y)
            self.assertAlmostEqual(s.norm(), 1.)

    def test_pdf(self):
        """Test the one-dimensional azimuthal response underlying pdf.
        """
        plt.figure('One-dimensional pdf')
        phi = numpy.linspace(-numpy.pi, numpy.pi, 100)
        for modulation in numpy.linspace(1, 0, 5):
            pdf = self.generator.pdf(phi, modulation)
            plt.plot(phi, pdf, label='$m = %.2f$' % modulation)
            spline = xInterpolatedUnivariateSplineLinear(phi, pdf)
            norm = spline.norm()
            self.assertTrue(abs(norm - 1.) < 1e-5,
                            'Normalization is %.3e' % norm)
        plt.axis([-numpy.pi, numpy.pi, None, None])
        plt.xlabel('$\\phi$ [rad]')
        plt.ylabel('pdf($\\phi$)')
        plt.legend(bbox_to_anchor=(0.88, 0.92))

    def test_cdf(self):
        """Test the one-dimensional azimuthal response underlying pdf.
        """
        plt.figure('One-dimensional cdf')
        phi = numpy.linspace(-numpy.pi, numpy.pi, 100)
        for modulation in numpy.linspace(1, 0, 5):
            cdf = self.generator.cdf(phi, modulation)
            plt.plot(phi, cdf, label='$m = %.2f$' % modulation)
            spline = xInterpolatedUnivariateSplineLinear(phi, cdf)
            self.assertTrue(abs(spline(phi[0]) - 0) < 1e-5, 'cdf(0) = %.3e' %\
                            spline(phi[0]))
            self.assertTrue(abs(spline(phi[-1]) - 1.) < 1e-5, 'cdf(1) = %.3e' %\
                            spline(phi[-1]))
        plt.axis([-numpy.pi, numpy.pi, None, None])
        plt.xlabel('$\\phi$ [rad]')
        plt.ylabel('cdf($\\phi$)')
        plt.legend(bbox_to_anchor=(0.4, 0.92))

    def test_ppf(self, num_points=100000):
        """Test the bilinear interpolation of the 2-dimensional ppf.
        """
        phi = numpy.random.uniform(-numpy.pi, numpy.pi, size=num_points)
        m = numpy.random.uniform(0., 1., size=num_points)
        q = self.generator.cdf(phi, m)
        delta = abs(self.generator.ppf(q, m) - phi) / phi
        std_err = numpy.std(delta)
        max_err = max(delta)
        idx = numpy.argmax(delta)
        logger.debug('Testing ppf spline interpolation on %d point(s)...' %\
                     num_points)
        logger.debug('Standard relative error: %f' % std_err)
        logger.debug('Maximum relative error: %f (at phi=%.3f, m=%.3f)' %\
                     (max_err, phi[idx], m[idx]))
        self.assertTrue(std_err < 1.e-3)
        self.assertTrue(max_err < 1.e-2)

    def test_rvs(self, num_events=1000000, m=0.25, phi0=0.25 * numpy.pi):
        """Test the random number generation.
        """
        mod = numpy.full(num_events, m)
        phi = self.generator.rvs_phi(mod, phi0)
        binning = numpy.linspace(-numpy.pi, numpy.pi, 100)
        hist = xHistogram1d(binning).fill(phi)
        model = xModulationCurveRad()
        fit_histogram(model, hist)
        if sys.flags.interactive:
            plt.figure()
            hist.plot()
            model.plot()
            model.stat_box()
        plt.xlabel('$\\phi$ [rad]')
        plt.axis([-numpy.pi, numpy.pi, 0, None])
        ahat, mhat, phihat = model.parameter_values()
        dahat, dmhat, dphihat = model.parameter_errors()
        mpull = (mhat - m) / dmhat
        phipull = (phihat - phi0) / dphihat
        self.assertTrue(abs(mpull) < 3.)
        self.assertTrue(abs(phipull) < 3.)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

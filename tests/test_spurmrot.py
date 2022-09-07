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

import numpy
from scipy.optimize import curve_fit

from ixpeobssim.core.hist import xHistogram1d, xHistogram2d
from ixpeobssim.irf.modf import xAzimuthalResponseGenerator
from ixpeobssim.evt.spurmrot import delta_phi_ampl, delta_phi_stokes,\
    correct_stokes_parameters, stokes_rotation_angle
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.core.modeling import xModulationCurveRad, xPowerLaw
from ixpeobssim.core.fitting import fit, fit_histogram


if sys.flags.interactive:
    plt.ion()


class TestSpurmrot(unittest.TestCase):

    """Unit test for the spurious modulation correction.
    """

    @classmethod
    def setUpClass(cls):
        """
        """
        cls.rnd = xAzimuthalResponseGenerator()

    @staticmethod
    def phi_to_stokes(phi, amplitude=2.):
        """Convenience function to convert a set of azimuthal angles (and, possibly,
        amplitudes) to Stokes parameters.
        """
        q = amplitude * numpy.cos(2. * phi)
        u = amplitude * numpy.sin(2. * phi)
        return q, u

    def generate_phi_sample(self, amplitude, phase, num_events):
        """Generate a sample of phi values, and fellow q and u Stokes parameters.
        """
        phi = self.rnd.rvs_phi(numpy.full(num_events, amplitude), phase)
        q, u = self.phi_to_stokes(phi)
        return phi, q, u

    def correct_phi(self, phi, q0, u0):
        """Run the default correction for a sample of phi values.
        """
        q, u = self.phi_to_stokes(phi)
        delta_phi = delta_phi_stokes(phi, q0, u0)
        phi_corr = phi + delta_phi
        q_corr, u_corr = self.phi_to_stokes(phi_corr)
        delta_q = q_corr - q
        delta_u = u_corr - u
        return delta_phi, phi_corr, q_corr, u_corr, delta_q, delta_u

    def _test_formulae_base(self, amplitude=0.1, phase=0., num_events=100000):
        """Test that the the three formulations for the correction
        are numerically equivalent.
        """
        phi, q, u = self.generate_phi_sample(amplitude, phase, num_events)
        q0, u0 = self.phi_to_stokes(phase, amplitude)
        delta_phi1 = delta_phi_ampl(phi, amplitude, phase)
        delta_phi2 = delta_phi_stokes(phi, q0, u0)
        # Compare the phi rotation in the amplitude/phase and Stokes parameter
        # flavors.
        self.assertTrue(numpy.allclose(delta_phi1, delta_phi2))
        phi_corr = phi + delta_phi1
        q_corr1, u_corr1 = self.phi_to_stokes(phi_corr)
        q_corr2, u_corr2 = correct_stokes_parameters(q, u, q0, u0)
        # Compare the corrected Stokes parameters from the corrected phi vs.
        # the direct rotation in Stokes space.
        self.assertTrue(numpy.allclose(q_corr1, q_corr2))
        self.assertTrue(numpy.allclose(u_corr1, u_corr2))

    def test_formulae(self):
        """Test the basic formulae for three different phase values.
        """
        self._test_formulae_base(0.1, 0.)
        self._test_formulae_base(0.1, numpy.radians(30.))
        self._test_formulae_base(0.1, numpy.radians(45.))
        self._test_formulae_base(0.1, numpy.radians(60.))
        self._test_formulae_base(0.1, numpy.radians(90.))

    def _test_sample_correction_base(self, amplitude=0.1, phase=0., num_events=1000000):
        """Base method to test the sample correction with no errors.
        """
        phi, q, u = self.generate_phi_sample(amplitude, phase, num_events)
        q0, u0 = self.phi_to_stokes(phase, amplitude)
        delta_phi, phi_corr, q_corr, u_corr, delta_q, delta_u = self.correct_phi(phi, q0, u0)
        plt.figure('Phi correction (amplitude=%.2f, phase=%.0f)' % (amplitude, numpy.degrees(phase)))
        phi_binning = numpy.linspace(-numpy.pi, numpy.pi, 250)
        h1 = xHistogram1d(phi_binning, xlabel=r'$\phi$ [$^\circ$]').fill(phi)
        h2 = xHistogram1d(phi_binning, xlabel=r'$\phi$ [$^\circ$]').fill(phi_corr)
        h1.plot()
        h2.plot()
        model = fit_histogram(xModulationCurveRad(), h2)
        m, dm = model.parameter_value('Modulation'), model.parameter_error('Modulation')
        self.assertTrue(abs(m) <= 5. * dm)
        model.plot()
        model.stat_box()

    @unittest.skip('Unnecessary for the moment...')
    def test_sample_correction(self):
        """Test the sample correction with no errors for three different phase values.
        """
        self._test_sample_correction_base(0.1, 0.)
        self._test_sample_correction_base(0.1, numpy.radians(30.))
        self._test_sample_correction_base(0.1, numpy.radians(45.))
        self._test_sample_correction_base(0.1, numpy.radians(60.))
        self._test_sample_correction_base(0.1, numpy.radians(90.))

    def test_average_stokes_rotation(self, amplitude=0.1, phase=numpy.radians(45.), num_events=1000000):
        """
        """
        phi, q, u = self.generate_phi_sample(amplitude, phase, num_events)
        q0, u0 = self.phi_to_stokes(phase, amplitude)
        delta = stokes_rotation_angle(q, u, q0, u0)
        delta2 = delta**2.
        print('(q0, u0) = (%.5f, %.5f)' % (q0, u0))
        print('Average (q, u) = (%.5f, %.5f)' % (q.mean(), u.mean()))
        print('Average (q^2, u^2) = (%.5f, %.5f)' % ((q**2.).mean(), (u**2.).mean()))
        print('Average qu = %.5f' % (q * u).mean())
        print('Average delta: %.5f' % delta.mean())
        print('Average delta^2: %.5f' % delta2.mean())
        print(0.5 * q0**2. + u0**2 * (1.5 - q0))
        plt.figure('Stokes rotation angle')
        binning = numpy.linspace(-3. * amplitude, 3. * amplitude, 250)
        h = xHistogram1d(binning).fill(delta)
        h.plot()
        setup_gca(xlabel='Stokes rotation angle [rad]')

    def test_errors(self, amplitude=0.1, phase=numpy.radians(0.), num_events=1000000):
        """
        """
        phi, q, u = self.generate_phi_sample(amplitude, phase, num_events)
        q0, u0 = self.phi_to_stokes(phase, amplitude)
        phi_binning = numpy.linspace(-numpy.pi, numpy.pi, 100)
        binning = numpy.linspace(-5. * amplitude, 5. * amplitude, 250)
        q_bias = []
        u_bias = []
        sigma_values = numpy.linspace(0., 0.1, 11)
        for sigma in sigma_values:
            qm = numpy.random.normal(q0, 0., size=phi.shape)
            um = numpy.random.normal(u0, sigma, size=phi.shape)
            delta_phi, phi_corr, q_corr, u_corr, delta_q, delta_u = self.correct_phi(phi, qm + sigma**2, um)
            plt.figure('Delta q')
            h = xHistogram1d(binning).fill(delta_q)
            bias = delta_q.mean() + q0
            q_bias.append(bias)
            h.plot(label=r'$\sigma = %.3f$, bias = %.3e' % (sigma, bias))
            setup_gca(xlabel=r'$\Delta q$', legend=True)
            plt.figure('Delta q 2d (sigma = %.3f)' % sigma)
            h = xHistogram2d(phi_binning, binning).fill(phi, delta_q)
            h.plot(logz=True)
            setup_gca(xlabel=r'$\phi$ [rad]', ylabel=r'$\Delta q$', grids=True)
            plt.figure('Delta u')
            h = xHistogram1d(binning).fill(delta_u)
            bias = delta_u.mean() + u0
            u_bias.append(bias)
            h.plot(label=r'$\sigma = %.3f$, bias = %.3e' % (sigma, bias))
            setup_gca(xlabel=r'$\Delta u$', legend=True)
        plt.figure('Bias')
        q_bias = numpy.array(q_bias)
        u_bias = numpy.array(u_bias)
        plt.plot(sigma_values, q_bias, 'o')

        def f(x, norm):
            return norm * x**2.

        popt, pcov = curve_fit(f, sigma_values, q_bias)
        x = numpy.linspace(0., sigma_values.max(), 100)
        plt.plot(x, f(x, *popt))
        print(q0, u0, popt)

        plt.plot(sigma_values, u_bias, 'o')

    def _test_errors_old(self, amplitude=0.1, phase=numpy.radians(0.), stokes_sigma=0.02, num_events=1000000):
        """
        """
        phi, q, u = self.generate_phi_sample(amplitude, phase, num_events)
        q0, u0 = self.phi_to_stokes(phase, amplitude)
        qgrid = numpy.linspace(q0 - 3. * stokes_sigma, q0 + 3. * stokes_sigma, 7)
        ugrid = numpy.linspace(u0 - 3. * stokes_sigma, u0 + 3. * stokes_sigma, 7)
        for _q in qgrid:
            for _u in ugrid:
                delta_phi, phi_corr, q_corr, u_corr, delta_q, delta_u = self.correct_phi(phi, _q, _u)
                print('(%.5f, %.5f) -> (%.6f, %.6f)' % (_q, _u, delta_q.mean() + _q, delta_u.mean() + _u))
            print()




if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

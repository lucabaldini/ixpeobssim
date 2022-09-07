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

from ixpeobssim.srcmodel.spectrum import power_law, cutoff_power_law
from ixpeobssim.srcmodel.spectrum import pl_integral, integral_flux
from ixpeobssim.srcmodel.spectrum import xSourceSpectrum, xCountSpectrum
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.irf import DEFAULT_IRF_NAME
from ixpeobssim.irf.caldb import irf_file_path
from ixpeobssim.irf.arf import xEffectiveArea
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()


class TestModels(unittest.TestCase):

    """Unit tests for the spectral models.
    """

    @staticmethod
    def _spectrum_spline(spectrum):
        """
        """
        fmt = dict(xlabel='Energy [keV]',
                   ylabel='dN/dE [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]')
        energy = numpy.linspace(1., 10., 200)
        return xInterpolatedUnivariateSpline(energy, spectrum(energy), **fmt)

    def test_power_law(self, pivot=5.):
        """
        """
        norm, index = 1., 2.
        pl = power_law(norm, index)
        self.assertAlmostEqual(norm * pivot**(-index), pl(pivot))
        s = self._spectrum_spline(pl)
        plt.figure('Power law')
        s.plot(logx=True, logy=True)

    def test_integral_flux(self):
        """
        """
        norm, index = 1., 2.
        emin, emax = 2., 8.
        spectrum = power_law(norm, index)
        val1 = pl_integral(norm, index, emin, emax)
        val2 = integral_flux(spectrum, emin, emax)
        self.assertAlmostEqual(val1, val2)



class TestSourceSpectrum(unittest.TestCase):

    """Unit tests for the xSourceSpectrum class.
    """

    def _basic_pl_test(self, column_density=0., redshift=0.):
        """Run a basic test with a power-law spectrum (where both the
        normalization and the index are time-dependent), and arbitary column
        density and redshift.
        """
        name = 'PL source spectrum (nH = %.2e, z = %.3f)' %\
               (column_density, redshift)
        # Create the xSourceSpectrum object.
        tmax = 10000.
        energy = numpy.linspace(1., 10., 181)
        time_ = numpy.linspace(0., tmax, 101)
        norm = lambda t: 1. + t / tmax
        index = lambda t: 2. - 0.5 * t / tmax
        pl = power_law(norm, index)
        spec = xSourceSpectrum(energy, time_, pl, column_density, redshift)
        # Basic two-dimentional plot of the thing.
        plt.figure(name)
        spec.plot(logz=True)
        # Time slices.
        plt.figure('%s time slices' % name)
        for _t in [0, 0.5 * tmax, tmax]:
            slice = spec.time_slice(_t)
            slice.plot(logx=True, logy=True, label='Slice @ t = %.1f s' % _t)
        plt.legend()
        # Energy slices.
        plt.figure('%s energy slices' % name)
        for _E in [1., 2., 8.]:
            slice = spec.energy_slice(_E)
            slice.plot(label='Slice @ E = %.1f keV' % _E)
        plt.legend()
        # Time-averaged spectrum.
        plt.figure('%s time-averaged spectrum' % name)
        spec.build_time_average().plot(logx=True, logy=True)
        # Light light_curve.
        plt.figure('%s light curve' % name)
        spec.build_light_curve().plot()

    def test_simple(self):
        """Test with a power-law, no interstellar absoprtion nor redhsift.
        """
        self._basic_pl_test(0., 0.)

    def test_absorbed(self):
        """
        """
        self._basic_pl_test(1.e22, 0.)



class TestCountSpectrum(unittest.TestCase):

    """Unit tests for the xCountSpectrum class.
    """

    @classmethod
    def setUpClass(cls, du_id=1):
        """Setup---here we essentially load the effective area.
        """
        arf_file_path = irf_file_path(DEFAULT_IRF_NAME, du_id, 'arf')
        cls.aeff = xEffectiveArea(arf_file_path)

    def test_naive(self):
        """Test a simple time-independent count spectrum.
        """
        norm = 1.
        index = 2.
        t = numpy.linspace(0., 1000., 101)
        pl = power_law(norm, index)
        cspec = xCountSpectrum(pl, self.aeff, t)
        plt.figure('PL count spectrum')
        cspec.plot()


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

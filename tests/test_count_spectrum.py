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
from ixpeobssim.irf import DEFAULT_IRF_NAME
from ixpeobssim.irf.caldb import irf_file_path
from ixpeobssim.irf.arf import xEffectiveArea
from ixpeobssim.srcmodel.spectrum import xCountSpectrum
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()

"""We explictely set the random seed to have reproducible results.
"""
numpy.random.seed(0)


class TestCountSpectrum(unittest.TestCase):

    """Unit test for xCountSpectrum.
    """

    @classmethod
    def setUpClass(cls, du_id=1):
        """Setup---here we essentially create the effective area.
        """
        arf_file_path = irf_file_path(DEFAULT_IRF_NAME, du_id, 'arf')
        cls.aeff = xEffectiveArea(arf_file_path)
        cls.interactive = sys.flags.interactive

    def test_power_law_stationary(self):
        """Test a time-independent power law.

        This creates a count spectrum with no time dependence, i.e., with
        only two (identical) interpolating time points on the auxiliary
        (time) axis. The power-law parameters (C and gamma) are constant
        and the underlying xUnivariateAuxGenerator looks like.

        .. image:: ../figures/test_power_law_stationary_2d.png

        Then a vertical slice (i.e., an interpolated linear spline) is taken
        in the middle of the auxiliary axis and the y-values of the spline are
        compared with the direct product of the effective area and the input spectrum.

        (Note that in general the two are technically not the same thing, as
        going from the count spectrum to the slice we do interpolate in time,
        although in this particular case the interpolation is trivial).
        If everything goes well, they should be on top of each other.
        The figure below is also showing the original power-law spectrum
        multiplied by the peak effective area.

        .. image:: ../figures/test_power_law_stationary_slice.png

        """
        tmin = 0.
        tmax = 100.
        tref = 0.5 * (tmin + tmax)
        C = 1.
        Gamma = 2.

        def powerlaw(E, t):
            """Function defining a time-dependent energy spectrum.
            """
            return C*numpy.power(E, -Gamma)

        plt.figure('Stationary PL count spectrum')
        _t = numpy.linspace(tmin, tmax, 5)
        count_spectrum = xCountSpectrum(powerlaw, self.aeff, _t)
        count_spectrum.plot()

        plt.figure('Stationary PL count spectrum slice @ %.3f s' % tref)
        ref_slice = count_spectrum.slice(tref)
        _x = self.aeff.x
        _y = C * numpy.power(_x, -Gamma) * self.aeff.y.max()
        plt.plot(_x, _y, '-', label='Original power-law spectrum')
        _y = C * numpy.power(_x, -Gamma) * self.aeff.y
        _mask = _y > 0.
        _x = _x[_mask]
        _y = _y[_mask]
        delta = abs((_y - ref_slice(_x)) / _y)
        delta_max = delta.max()
        plt.plot(_x, _y, 'o', label='Direct convolution with aeff')
        ref_slice.plot(logx=True, logy=True, label='xCountSpectrum output')
        plt.text(0.1, 0.1, 'Max. difference = %.3e' % delta_max,
                 transform=plt.gca().transAxes)
        plt.legend(bbox_to_anchor=(0.75, 0.5))
        plt.axis([0.6, 10, None, None])
        self.assertTrue(delta_max < 1e-3, 'max deviation %.9f' % delta_max)

        plt.figure('Time-integrated spectrum')
        spec = count_spectrum.build_time_integral()
        spec.plot()

    def test_power_law_variable(self):
        """Test a time-dependent power law.

        This creates a time-dependent count spectrum, where the two parameters
        of the underlying power law (C and Gamma) vary linearly with time, in
        opposite direction, between 1 and 2.

        (Beware: this does not mean that you can interpolate linearly between
        the two time extremes, as both parameters vary at the same time and
        the spectral shape does not evolve linearly with time---we're sampling
        the time axis with 100 points).
        The underlying xUnivariateAuxGenerator looks like.

        .. image:: ../figures/test_power_law_variable_2d.png

        Then a vertical slice (i.e., an interpolated linear spline) is taken
        in the middle of the auxiliary axis and the y-values of the spline
        are compared with the direct product of the effective area and the
        count spectrum (evaluated at the same time). If everything goes well,
        they should be on top of each other. The figure below is also showing
        the orignal power-law spectrum multiplied by the peak effective area.

        .. image:: ../figures/test_power_law_variable_slice.png

        Finally, we do test the light-curve building by comparing it with the
        values from a direct intergration of the vertical slices on a
        fixed-spacing grid. Note that, since the normalization increases with
        time and the spectral index becomes harder, the light-curve increases
        more than linearly.

        .. image:: ../figures/test_power_law_variable_lc.png

        """
        tmin = 0.
        tmax = 100.
        tref = 0.5*(tmin + tmax)

        def C(t):
            """Time-dependent C---equals to 1 @ tmin and 2 @ tmax.
            """
            return 1. + (t - tmin) / (tmax - tmin)

        def Gamma(t):
            """Time-dependent C---equals to 2 @ tmin and 1 @ tmax.
            """
            return 2. - (t - tmin) / (tmax - tmin)

        def powerlaw(E, t):
            """Function defining a time-dependent energy spectrum.
            """
            return C(t) * numpy.power(E, -Gamma(t))

        plt.figure('Variable PL count spectrum')
        _t = numpy.linspace(tmin, tmax, 100)
        count_spectrum = xCountSpectrum(powerlaw, self.aeff, _t)
        count_spectrum.plot()

        plt.figure('Variable PL count spectrum slice @ %.3f s' % tref)
        ref_slice = count_spectrum.slice(tref)
        _x = self.aeff.x
        _y = self.aeff.y.max()*C(tref)*numpy.power(_x, -Gamma(tref))
        plt.plot(_x, _y, '-', label='Original power-law spectrum')
        _y = C(tref)*numpy.power(_x, -Gamma(tref))*self.aeff.y
        _mask = _y > 0.
        _x = _x[_mask]
        _y = _y[_mask]
        delta = abs((_y - ref_slice(_x))/_y).max()
        self.assertTrue(delta < 1e-3, 'max deviation %.9f' % delta)
        plt.plot(_x, _y, 'o', label='Direct convolution with aeff')
        ref_slice.plot(logx=True, logy=True, label='xCountSpectrum output')
        plt.text(0.1, 0.1, 'Max. difference = %.3e' % delta,
                 transform=plt.gca().transAxes)
        plt.legend(bbox_to_anchor=(0.75, 0.5))
        plt.axis([0.6, 10, None, None])

        plt.figure('Variable PL light curve')
        _x = numpy.linspace(tmin, tmax, 33)
        _y = []
        for _xp in _x:
            _y.append(count_spectrum.slice(_xp).norm())
        _y = numpy.array(_y)
        plt.plot(_x, _y, 'o', label='Direct integral flux values')
        delta = abs((_y - count_spectrum.light_curve(_x))/_y).max()
        self.assertTrue(delta < 1e-3, 'max deviation %.9f' % delta)
        count_spectrum.light_curve.plot(label='xCountSpectrum light-curve')
        plt.legend(bbox_to_anchor=(0.65, 0.75))
        plt.text(0.5, 0.1, 'Max. difference = %.3e' % delta,
                 transform=plt.gca().transAxes)

    def test_power_law_rvs(self, num_events=1000000):
        """Test the generation of event energies from a count power-law
        spectrum convoluted with the effective area.

        This turned out to be more tricky than we anticipated. Since the
        convolution of the source spectrum with the effective area falls
        pretty quickly at high energy, there's typically very few events above
        a few keV and, as a consequence, the slope of the corresponding ppf is
        fairly steep close to one.

        .. image:: ../figures/test_power_law_rvs_ppf.png

        This implies that the ppf in each slice must be properly sampled
        close to 1 (initial tests showed, e.g., that a uniform grid with
        100 points between 0 and 1 was not enough to throw meaningful random
        numbers for a typical power-law source spectrum). This is particularly
        true for soft spectral indices---which is why we picked `Gamma = 3.`
        for this test.

        .. image:: ../figures/test_power_law_rvs_counts.png

        """
        tmin = 0.
        tmax = 100.
        tref = 0.5 * (tmin + tmax)
        C = 1.
        Gamma = 3.

        def powerlaw(E, t):
            """Function defining a time-dependent energy spectrum.
            """
            return C*numpy.power(E, -Gamma)

        plt.figure('Count spectrum ppf slice @ %.3f s' % tref)
        _t = numpy.linspace(tmin, tmax, 5)
        count_spectrum = xCountSpectrum(powerlaw, self.aeff, _t)
        count_spectrum.ppf.hslice(tref).plot(overlay=True)

        plt.figure('Count spectrum random variates')
        ref_slice = count_spectrum.slice(tref)
        _time = numpy.full(num_events, tref)
        _energy = count_spectrum.rvs(_time)
        _binning = numpy.linspace(self.aeff.xmin(), self.aeff.xmax(), 100)
        obs, bins, patches = plt.hist(_energy, bins=_binning,
                                      histtype='step', label='Random energies')
        plt.yscale('log')
        # We want to overlay the reference count-spectrum slice, normalized
        # to the total number of events simulated.
        bin_width = (bins[1] - bins[0])
        scale = num_events * bin_width / ref_slice.norm()
        _x = 0.5*(bins[:-1] + bins[1:])
        _y = scale*ref_slice(_x)
        plt.plot(_x, _y, label='Underlying pdf')
        # And, for the chisquare, we do correctly integrate the slice in each
        # energy bin, rather than evaluating it at the bin center.
        exp = []
        scale = num_events/ref_slice.norm()
        for _emin, _emax in zip(bins[:-1], bins[1:]):
            exp.append(scale*ref_slice.integral(_emin, _emax))
        exp = numpy.array(exp)
        _mask = exp > 0.
        exp = exp[_mask]
        obs = obs[_mask]
        chi2 = ((exp - obs)**2/exp).sum()
        ndof = len(obs)
        chi2_min = ndof - 3*numpy.sqrt(2*ndof)
        chi2_max = ndof + 3*numpy.sqrt(2*ndof)
        self.assertTrue(chi2 > chi2_min, 'chisquare too low (%.2f/%d)' %\
                        (chi2, ndof))
        self.assertTrue(chi2 < chi2_max, 'chisquare too high (%.2f/%d)' %\
                        (chi2, ndof))
        plt.text(0.5, 0.1, '$\chi^2$/ndof = %.2f/%d' % (chi2, ndof),
                 transform=plt.gca().transAxes)
        plt.legend(bbox_to_anchor=(0.85, 0.75))



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

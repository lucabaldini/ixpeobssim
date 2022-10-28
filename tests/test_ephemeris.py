#!/usr/bin/env python
#
# Copyright (C) 2020, the ixpeobssim team.
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
import datetime
import sys
import os

import numpy
from matplotlib.dates import date2num
from mpl_toolkits.mplot3d import Axes3D

from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.core.rand import xUnivariateGenerator
from ixpeobssim.srcmodel.ephemeris import xEphemeris, xOrbitalEphemeris
from ixpeobssim.srcmodel.ephemeris import time_list, read_par_file
from ixpeobssim.srcmodel.ephemeris import get_eccentric_anomaly, get_earth_barycentric_ephemeris
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()


class TestEphemeris(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Crab ephemeris:
        * nu0 (Hz)                : 29.8003951530036
        * nudot0 (10^-10 Hz s^-1) : -3.73414
        * nuddot (10^-20 Hz s^-2) : 1.18
        """
        cls.crab_ephemeris = xEphemeris(0., 29.8003951530036, -3.73414e-10, 1.18e-20)
        print(cls.crab_ephemeris)

    def test_basic(self, met0=0., nu0=2.):
        """Basic test
        """
        ephem = xEphemeris(met0, nu0)
        print(ephem)
        print(ephem.dict())
        self.assertAlmostEqual(ephem.met_to_phase(0.), 0.)
        self.assertAlmostEqual(ephem.met_to_phase(2. * numpy.pi / nu0), 2. * numpy.pi)
        met = numpy.linspace(met0, met0 + 1000., 100)
        phase = ephem.met_to_phase(met)
        self.assertTrue(numpy.allclose(phase, met * nu0))

    def test_spline(self, met0=10., nu0=0.033, duration=1.e6):
        """
        """
        ephem = xEphemeris(met0, nu0)
        plt.figure('Ephemeris spline')
        ephem.phase_spline(met0, duration).plot()

    def test_spline_inverse(self, met0=10., nu0=0.033, duration=1.e6):
        """
        """
        ephem = xEphemeris(met0, nu0)
        plt.figure('Ephemeris spline inverse')
        ephem.phase_spline_inverse(met0, duration).plot()

    def test_roundtrip(self, duration=1.e7, padding=1.e10):
        """
        """
        met0 = self.crab_ephemeris.met0
        for start_met in (met0, met0 + padding, met0 - padding):
            met1 = numpy.linspace(start_met, start_met + duration, 100000)
            phase = self.crab_ephemeris.met_to_phase(met1)
            met2 = self.crab_ephemeris.phase_spline_inverse(start_met, duration)(phase)
            self.assertTrue(numpy.allclose(met2, met1))
            met3 = self.crab_ephemeris.phase_to_met(phase, start_met, duration)
            self.assertTrue(numpy.allclose(met3, met1))

    def _test_ephemeris_rvs(self, ephem, label, start_met=1.e6, duration=10000,
                            num_events=1000000, test_chi2=True, atol=1.e-8):
        """
        """
        num_bins = 250
        x = numpy.linspace(0., 1., num_bins)
        y = numpy.sin(x * 2 * numpy.pi)**2.
        pulse_profile = xUnivariateGenerator(x, y)
        pulse_phase, met = ephem.rvs(pulse_profile, start_met, duration, num_events)
        self.assertEqual(len(pulse_phase), num_events)
        self.assertEqual(len(met), num_events)
        self.assertTrue(met.min() >= start_met)
        self.assertTrue(met.max() <= start_met + duration)
        # Make sure that the pulse phase and the MET are consistent with each
        # other, and that the ephemeris conversion functions do roundtrip.
        # The error message is meant to help debugging rounding errors.
        _pp = ephem.fold(met, start_met)
        args = pulse_phase, _pp, _pp - pulse_phase
        msg = 'Roundtrip error...\nPulse phase = %s\nFolded MET = %s\nDifference = %s' % args
        self.assertTrue(numpy.allclose(pulse_phase, _pp, atol=atol), msg)
        # And now on with the plots.
        plt.figure('Test %s phase' % label)
        binning = numpy.linspace(0., 1., num_bins + 1)
        hist = xHistogram1d(binning, xlabel='Pulse phase').fill(pulse_phase)
        hist.plot()
        if test_chi2:
            # Calculate a chisquare of the pulse profile and make sure that it
            # makes sense.
            scale_factor = num_events * 2. / num_bins
            x = hist.bin_centers(0)
            y = hist.content
            mask = y > 0.
            x = x[mask]
            y = y[mask]
            chisq = ((y - scale_factor * pulse_profile(x))**2. / y).sum()
            delta = abs(chisq - num_bins) / numpy.sqrt(num_bins)
            self.assertTrue(delta < 10.)
            pulse_profile.plot(scale=scale_factor)
        plt.figure('Test %s met' % label)
        binning = numpy.linspace(start_met, start_met + duration, num_bins + 1)
        hist = xHistogram1d(binning, xlabel='MET [s]').fill(met)
        hist.plot()

    def test_ephermeris_rvs1(self):
        """Simple test, no frequency derivatives.
        """
        ephem = xEphemeris(-1.e6, 0.13)
        self._test_ephemeris_rvs(ephem, 'ephem 1')

    def test_ephermeris_rvs2(self):
        """Standard test.
        """
        ephem = xEphemeris(0., 0.13, 1.e-5, 1.e-9)
        self._test_ephemeris_rvs(ephem, 'ephem 2', start_met=10.)

    def test_ephermeris_rvs3(self):
        """Specific test for the fractional part of the last period.
        """
        ephem = xEphemeris(0., 1.)
        self._test_ephemeris_rvs(ephem, 'ephem 3', start_met=10., duration=1.5, test_chi2=False)

    def test_ephermeris_rvs4(self):
        """Test with the Crab ephemeris, and the start met ~20 years after the
        ephermeris met0.
        """
        self._test_ephemeris_rvs(self.crab_ephemeris, 'ephem Crab', start_met=1.e9,
                                 duration=1.e5, atol=1.e-4)

    def test_time_list(self, num_events=1000000, duration=1000000.):
        """Not sure what this is really testing...
        """
        file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'parfiles', 'SAXJ1808.4-3658.par')
        ephem = xOrbitalEphemeris.from_file(file_path)
        start_met = ephem.met0
        phase = numpy.linspace(0., 1., num_events)
        inv_spline = ephem.phase_spline_inverse(start_met, duration)
        num_periods = int(duration * ephem.nu0)
        p = numpy.random.randint(0, num_periods, num_events)
        met = inv_spline(p + phase)
        times = time_list(inv_spline, start_met, ephem, num_events, duration)

    @unittest.skip('Tempoarily disable due to a weird failure in https://github.com/lucabaldini/ixpeobssim/pull/642')
    def test_get_earth_barycentric_ephemeris(self):
        """
        """
        met = numpy.linspace(0, 1e6, 1000)
        ebe = get_earth_barycentric_ephemeris(met)
        fig = plt.figure('Earth barycentric ephemeris')
        ax = fig.gca(projection='3d')
        ax.plot(ebe[0], ebe[1], ebe[2])

    def test_get_eccentric_anomaly(self):
        """
        """
        met = numpy.linspace(0, 1e6, 1000)
        ephem = xEphemeris(10., 0.033)
        ephem = xOrbitalEphemeris(ephem.met0, ephem.nu0, 10., 100., 0.)
        periods = 1000
        samples = 100
        ea = get_eccentric_anomaly(met, ephem, periods, samples, plot=True)

    def test_read_par_file(self):
        """Read a .par file.
        """
        file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'parfiles', 'SAXJ1808.4-3658.par')
        params = read_par_file(file_path)
        print(params)

    def test_orb_ephem(self):
        """
        """
        file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'parfiles', 'SAXJ1808.4-3658.par')
        ephem = xOrbitalEphemeris.from_file(file_path)
        print(ephem)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

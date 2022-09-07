#!/urs/bin/env python
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

import sys
import unittest

import numpy

from ixpeobssim.instrument.traj import xIXPETrajectory, xSAABoundary, xObservationTimeline
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.time_ import met_to_mjd, seconds_to_days, years_to_seconds
from ixpeobssim.utils.profile import xChrono
from ixpeobssim.utils.logging_ import logger

if sys.flags.interactive:
    plt.ion()



class TestTrajectory(unittest.TestCase):

    """Unit test for xIXPETrajectory.
    """

    @classmethod
    def setUpClass(cls):
        """Setup the test.
        """
        cls.start_met = 2200.
        cls.duration = 25800.
        cls.stop_met = cls.start_met + cls.duration
        cls.step = 10.
        cls.met = numpy.arange(cls.start_met, cls.stop_met, cls.step)
        cls.trajectory = xIXPETrajectory()
        cls.saa = xSAABoundary()
        print(cls.saa)

    def test_trajectory(self):
        """Basic orbit plot.
        """
        plt.figure('Satellite trajectory')
        lon, lat = self.trajectory.position(self.met)
        # Split the trajectory into segements to avoid annoying horizontal lines
        # in the plot when passing from 180 to -180 in longitude.
        idx = numpy.where(numpy.diff(lon) < 0.)[0] + 1
        lon = numpy.split(lon, idx)
        lat = numpy.split(lat, idx)
        for _lon, _lat in zip(lon, lat):
            plt.plot(_lon, _lat, color='black')
            mask = self.saa.contains(_lon, _lat)
            plt.plot(_lon[mask], _lat[mask], color='lightgray')
        self.saa.plot()
        setup_gca(xmin=-180., xmax=180., ymin=-6., ymax=6., grids=True,
                  xlabel='Longitude [deg]', ylabel='Latitude [deg]')

    def test_saa(self, precision=0.001):
        """Calculation of the SAA epochs.
        """
        plt.figure('SAA')
        plt.plot(self.met, self.trajectory.in_saa(self.met))
        setup_gca(ymin=0., ymax=1.2, xlabel='MET [s]', ylabel='SAA flag',
                  xmin=self.start_met, xmax= self.stop_met)
        args = self.start_met, self.stop_met
        epochs = self.trajectory.saa_epochs(*args, precision=precision)
        print('SAA epochs: %s' % epochs)
        _test_function = lambda met: self.trajectory.in_saa(met)
        delta = precision
        for (t1, t2) in epochs:
            # This is the actual test: if we move by a quantity >= precision from
            # the epoch bounds we know exactly what we should expect.
            if t1 > self.start_met:
                self.assertFalse(_test_function(t1 - delta))
            self.assertTrue(_test_function(t1 + delta))
            self.assertTrue(_test_function(t2 - delta))
            if t2 < self.stop_met:
                self.assertFalse(_test_function(t2 + delta))
            plt.axvline(t1, ls='dashed', color='red')
            plt.axvline(t2, ls='dashed', color='blue')

    def test_saa_duration(self, duration=1000000):
        """Plot the distribution of the duration of the SAA passages.
        """
        logger.info('Testing SAA epoch duration over %s s...', duration)
        epochs = self.trajectory.saa_epochs(0., duration, 100)
        met = numpy.array(tuple(_met for epoch in epochs for _met in epoch))
        dist = met[2:-1:2] - met[0:-2:2]
        logger.info('Average distance between SAA entries: %.3f +- %.3f s',\
            dist.mean(), dist.std(ddof=1))
        delta = numpy.array([t2 - t1 for (t1, t2) in epochs])
        min_ = delta.min()
        max_ = delta.max()
        mean = delta.mean()
        rms = delta.std(ddof=1)
        logger.info('%d epoch(s) found, spanning %.3f---%.3f s.', len(epochs), min_, max_)
        plt.figure('SAA epoch duration')
        plt.hist(delta, 100)
        setup_gca(xlabel='SAA epoch duration [s]', ylabel='Entries per bin')

    def test_earth_occultation(self, ra=83.6330833, dec=22.0145000, precision=0.001):
        """Calculation of the Earth occultation epochs.
        """
        plt.figure('Earth occultation')
        plt.plot(self.met, self.trajectory.target_occulted(self.met, ra, dec))
        setup_gca(ymin=0., ymax=1.2, xlabel='MET [s]', ylabel='Earth occultation',
                  xmin=self.start_met, xmax= self.stop_met)
        args = self.start_met, self.stop_met, ra, dec
        epochs = self.trajectory.earth_occultation_epochs(*args, precision=precision)
        print('Earth occultation epochs: %s' % epochs)
        _test_function = lambda met: self.trajectory.target_occulted(met, ra, dec)
        delta = 2. * precision
        for (t1, t2) in epochs:
            # This is the actual test: if we move by a quantity >= precision from
            # the epoch bounds we know exactly what we should expect.
            if t1 > self.start_met:
                self.assertFalse(_test_function(t1 - delta))
            self.assertTrue(_test_function(t1 + delta))
            self.assertTrue(_test_function(t2 - delta))
            if t2 < self.stop_met:
                self.assertFalse(_test_function(t2 + delta))
            plt.axvline(t1, ls='dashed', color='red')
            plt.axvline(t2, ls='dashed', color='blue')

    def _gti_list(self, *args, **kwargs):
        """
        """
        gtis, _, _ = self.trajectory.gti_list(self.start_met, self.stop_met, *args, **kwargs)
        return gtis

    def test_gtis(self, ra=83.6330833, dec=22.0145000):
        """
        """
        gtis = self._gti_list(ra, dec, saa=False, occult=False)
        print(gtis)
        self.assertEqual(len(gtis), 1)
        start, stop = gtis[0]
        self.assertAlmostEqual(start, self.start_met)
        self.assertAlmostEqual(stop, self.stop_met)
        gtis = self._gti_list(ra, dec, saa=False)
        print(gtis)
        gtis = self._gti_list(ra, dec, occult=False)
        print(gtis)
        plt.figure('GTI')
        plt.plot(self.met, self.trajectory.in_saa(self.met), label='SAA')
        plt.plot(self.met, self.trajectory.target_occulted(self.met, ra, dec), label='Earth occultation')
        setup_gca(ymin=0., ymax=1.4, xlabel='MET [s]', xmin=self.start_met,
                  xmax= self.stop_met, legend=True)
        gtis = self._gti_list(ra, dec)
        print(gtis)
        y0 = 1.1
        for (t1, t2) in gtis:
            plt.hlines(y0, t1, t2, color='black')
            plt.plot([t1, t2], [y0, y0], 'o', color='black')

    def test_epochs(self):
        """
        """
        start_met = 0.
        stop_met = 1.
        epochs = []
        complement = self.trajectory.complement_epochs(epochs, start_met, stop_met)
        self.assertEqual(complement, [(start_met, stop_met)])
        complement = self.trajectory.complement_epochs(complement, start_met, stop_met)
        self.assertEqual(complement, [])

    def test_timeline(self, start_met=0., stop_met=100000., ra=45., dec=45.):
        """
        """
        timeline = xObservationTimeline(start_met, stop_met, ra, dec, True, True)
        timeline_gti_list = timeline.gti_list()
        print(timeline_gti_list)
        traj_gti_list, _, _ = timeline.trajectory.gti_list(start_met, stop_met, ra, dec)
        for gti1, gti2 in zip(timeline_gti_list, traj_gti_list):
            self.assertAlmostEqual(gti1, gti2)

    def test_scdata(self, start_met=1., stop_met=9999., ra=45., dec=45., time_step=30.):
        """
        """
        timeline = xObservationTimeline(start_met, stop_met, ra, dec, True, True)
        sc_data = timeline.sc_data(time_step)
        plt.figure('Observation ground track')
        plt.plot(sc_data['LON_GEO'], sc_data['LAT_GEO'], '.')
        setup_gca(xlabel='Longitude [deg]', ylabel='Latitude [deg]', grids=True)
        plt.figure('Observation flags')
        plt.plot(sc_data['MET'], sc_data['IN_SAA'], label='SAA')
        plt.plot(sc_data['MET'], sc_data['TARGET_OCCULT'], label='Target occulted')
        setup_gca(xlabel='MET [s]', ylabel='Flags', legend=True, grids=True,
                  xmin=start_met, xmax=stop_met, ymax=1.25)


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

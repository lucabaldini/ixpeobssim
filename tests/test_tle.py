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

import os

import sys
import unittest

import numpy
import skyfield.api

from ixpeobssim import IXPEOBSSIM_INSTRUMENT_DATA
from ixpeobssim.instrument.traj import xSkyfieldLoader, xIXPETrajectory, xTLE
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.time_ import met_to_jd

if sys.flags.interactive:
    plt.ion()


class TestTLE(unittest.TestCase):

    """Unit test for the TLE facilities.
    """

    @classmethod
    def setUpClass(cls, tle_file_name='tle_data_science.txt'):
        """Setup the test.
        """
        tle_file_path = os.path.join(IXPEOBSSIM_INSTRUMENT_DATA, tle_file_name)
        logger.info('Loading IXPE TLE data from %s...', tle_file_path)
        cls.satellites = {s.name: s for s in skyfield.api.load.tle_file(tle_file_path)}
        cls.timescale = xSkyfieldLoader().timescale(builtin=True)
        line1, line2 = xTLE.lines()
        ixpe = skyfield.api.EarthSatellite(line1, line2, 'IXPE', cls.timescale)
        cls.satellites['IXPE'] = ixpe
        print(xTLE._epoch())

    def _test_satellite(self, name, duratio=10000):
        """
        """
        satellite = self.satellites[name]
        print(satellite)
        met = numpy.linspace(0., 10000., 100)
        t = self.timescale.ut1_jd(met_to_jd(met))
        geocentric = satellite.at(t)
        subpoint = geocentric.subpoint()
        alt = subpoint.elevation.km
        lon = subpoint.longitude.degrees
        lat = subpoint.latitude.degrees
        plt.figure('%s altitude' % name)
        plt.plot(met, alt)
        setup_gca(xlabel='MET [s]', ylabel='Altitude [km]', grids=True)
        plt.figure('%s orbit' % name)
        plt.plot(lon, lat)
        setup_gca(xlabel='Longitude [degrees]', ylabel='Latitude [degrees]', grids=True)

    def test_agile(self, duration=10000):
        """
        """
        self._test_satellite('AGILE')

    def test_nustar(self):
        """
        """
        self._test_satellite('NUSTAR')

    def test_ixpe(self):
        """
        """
        self._test_satellite('IXPE')

    def test_ixpe_tle(self):
        """The IXPE TLE should look like
        1 49954U 21121A   21351.00640149  .00001120  00000-0  35770-4 0  9994
        2 49954   0.2300 281.7657 0011347 134.4260 303.9164 14.90740926  1166

        This TLE is coming from
        https://www.n2yo.com/satellite/?s=49954
        and was taken on December, 17, 2021.
        """
        line1, line2 = xTLE.lines()
        print(line1)
        print(line2)
        self.assertEqual(line1, '1 49954U 21121A   21351.00640149  .00001120  00000-0  35770-4 0  9994')
        self.assertEqual(line2, '2 49954   0.2300 281.7657 0011347 134.4260 303.9164 14.90740926  1166')




if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

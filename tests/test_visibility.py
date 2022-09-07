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
import os

import numpy
import matplotlib.dates

from ixpeobssim.instrument.traj import xIXPETrajectory
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.time_ import string_to_met_utc, met_to_num
from ixpeobssim.utils.logging_ import logger
from ixpeobssim import IXPEOBSSIM_TEST_DATA

if sys.flags.interactive:
    plt.ion()



class TestVisibility(unittest.TestCase):

    """Unit test for the source visibility.
    """

    @classmethod
    def setUpClass(cls):
        """Setup the test.
        """
        cls.trajectory = xIXPETrajectory()

    def test_sun_position(self):
        """From https://www.calsky.com/cs.cgi?cha=5&sec=1 set to April, 5, 2020

        Topocentric:
	    Coordinates including refraction
	    Altitude:    -36.852°             Azimuth:    335.157°    Direction: North-Northwest NNW
	    Apparent:    R.A.  0h 57m 26s   Dec. +   6° 08' 03"
        """
        met = string_to_met_utc('2020-04-05', lazy=True)
        ut = self.trajectory.met_to_ut(met)
        pos = self.trajectory.planets['earth'].at(ut).observe(self.trajectory.planets['sun'])
        print(pos.radec())

    def test_chandra_crab_visibility(self):
        """From https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/viewing/viewing.pl

        Viewing Results
        Input equatorial coordinates:
        Crab, resolved by SIMBAD (local cache) to
        [ 83.6331°, 22.0145° ], equinox J2000.0

        IXPE (planning)

        This mission is in the planning stage.

        *** VIEWING Version 3.4      run on 2020 May 13 ***
        for the period 2020 May 13 to 2022 May 14

        With IXPE (Sun angle range =  65-115):
        Observable between 2020 Aug 22           and 2020 Oct 12
        Observable between 2021 Feb 17           and 2021 Apr 09
        Observable between 2021 Aug 22           and 2021 Oct 12
        Observable between 2022 Feb 18           and 2022 Apr 09
        """
        ra, dec = 83.63308333, 22.0145
        file_path = os.path.join(IXPEOBSSIM_TEST_DATA, 'chandra_crab_visibility.txt')
        year, month, day, mjd, vis, pitch = numpy.loadtxt(file_path, unpack=True)
        year = year.astype(int)
        month = month.astype(int)
        day = day.astype(int)
        met = numpy.array([string_to_met_utc('%d-%d-%d' % (y, m, d), lazy=True)\
            for y, m, d in zip(year, month, day)])
        diff = self.trajectory.target_sun_angle(met, ra, dec) - pitch
        plt.figure('Pitch angle difference')
        plt.plot(met_to_num(met), diff)
        self.assertTrue(abs(diff).max() < 2.)
        formatter = matplotlib.dates.DateFormatter('%Y-%m-%d')
        plt.gca().xaxis.set_major_formatter(formatter)
        setup_gca(xlabel='Date and time', ylabel='Sun angle difference [degrees]',
                  grids=True)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

#!/usr/bin/env python
#
# Copyright (C) 2018, the ixpeobssim team.
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
import astropy.io.fits as fits

from ixpeobssim.utils.logging_ import logger
import ixpeobssim.config.toy_point_source as toy_point_source
from ixpeobssim.utils.time_ import string_to_met_utc
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.utils.environment import PYXSPEC_INSTALLED
if PYXSPEC_INSTALLED:
    import ixpeobssim.evt.xspec_ as xspec_



class TestObssim(unittest.TestCase):

    """Unit test for simulation and analysis pipeline(s).
    """

    def test_empty_event_list(self):
        """This is to verify that xpobssim does not crash with empty event lists.
        """
        pipeline.reset('toy_point_source', overwrite=True)
        evt_file_list = pipeline.xpobssim(duration=1.e-3)

    def test_toy_point_source(self):
        """Simulate a simple source with a power-law spectrum, fit the counts
        with XSPEC and make sure that the spectral parameters are consistent
        with the input setup.
        """
        pipeline.reset('toy_point_source', overwrite=True)
        evt_file_list = pipeline.xpobssim(duration=1000., deadtime=0.,
                                          saa=False, occult=False)
        pha_file_list = pipeline.xpbin(*evt_file_list, algorithm='PHA1')
        if not PYXSPEC_INSTALLED:
            return
        fit_output = pipeline.xpxspec(*pha_file_list, model='powerlaw', plot=False)
        target = (toy_point_source.pl_index, toy_point_source.pl_norm)
        self.assertTrue(xspec_.compare_fit_data(fit_output, target) == 0)

    @unittest.skip("Disengaged for the public release")
    def test_crab_pulsar(self):
        """Run a simple simulation with the Crab pulsar.
        """
        pipeline.reset('crab_pulsar', overwrite=True)
        kwargs = dict(startdate='2022-04-21', duration=1000., deadtime=0.,
                      saa=False, occult=False)
        evt_file_list = pipeline.xpobssim(**kwargs)
        pha_file_list = pipeline.xpbin(*evt_file_list, algorithm='PHA1')
        # Test pulsar timing
        hdu = fits.open(evt_file_list[0])
        t = hdu['EVENTS'].data.field('TIME')
        tmin = t.min()
        tmax = t.max()
        # Verify that the minimum event time is larger than start MET
        start_met = string_to_met_utc(kwargs.get('startdate'), lazy=True)
        self.assertTrue(tmin >= start_met)
        # Verify that the time difference between the last and the first event
        # is smaller than the duration
        self.assertTrue(tmax - tmin <= kwargs.get('duration'))
        # Verify that the time column is sorted.
        self.assertTrue((numpy.diff(t) >= 0).all())



if __name__ == '__main__':
    unittest.main()

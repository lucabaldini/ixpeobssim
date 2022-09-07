#!/urs/bin/env python
#
# Copyright (C) 2020--2022, the ixpeobssim team.
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

from ixpeobssim.binning.misc import xBinnedMap
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.core.hist import xHistogram2d
from ixpeobssim.utils.logging_ import logger
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
import ixpeobssim.config.toy_multiple_sources as srcmod

if sys.flags.interactive:
    plt.ion()


class TestPolarizationMapCube(unittest.TestCase):

    """Unit test for the binned polarization cubes.
    """

    @classmethod
    def setUpClass(cls):
        """Run a small simulation and create the polarization cubes.
        """
        pipeline.reset('toy_multiple_sources', overwrite=True)
        cls.file_list = pipeline.xpobssim(duration=10000., seed=13)

    def test_file(self):
        """
        """
        file_path = self.file_list[0]
        event_file = xEventFile(file_path)
        ra, dec = event_file.sky_position_data()
        mc_ra, mc_dec = event_file.sky_position_data(True)
        ra0, dec0 = numpy.mean(ra), numpy.mean(dec)
        mc_ra0, mc_dec0 = numpy.mean(mc_ra), numpy.mean(mc_dec)
        logger.info('Average measured sky position: (%.5f, %.5f)', ra0, dec0)
        logger.info('Average MC sky position: (%.5f, %.5f)', mc_ra0, mc_dec0)
        delta = 0.075
        num_bins = 100
        binning = (numpy.linspace(srcmod.ra - delta, srcmod.ra + delta, num_bins),
                   numpy.linspace(srcmod.dec - delta, srcmod.dec + delta, num_bins))
        hist = xHistogram2d(*binning).fill(ra, dec)
        plt.figure('Raw counts')
        hist.plot()
        setup_gca(grids=True)

    def test_data(self):
        """
        """
        file_list = pipeline.xpbin(*self.file_list, algorithm='CMAP')
        cmap = xBinnedMap.from_file_list(file_list)
        plt.figure('Measured map')
        cmap.plot()

    def test_mc(self):
        """
        """
        file_list = pipeline.xpbin(*self.file_list, algorithm='CMAP', mc=True, suffix='cmap_mc')
        cmap = xBinnedMap.from_file_list(file_list)
        plt.figure('MC map')
        cmap.plot()



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

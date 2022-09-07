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

"""Unit tests for srcmodel.img
"""


from __future__ import print_function, division

import sys
import os
import unittest

import numpy

from ixpeobssim import IXPEOBSSIM_CONFIG_FITS, IXPEOBSSIM_TEST_DATA
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.srcmodel.img import xFITSImage
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()



class TestImg(unittest.TestCase):

    """.
    """

    def test_toy_point_source(self):
        """Test with a toy point source.
        """
        file_path = os.path.join(IXPEOBSSIM_TEST_DATA, 'point_count_map.fits')
        image = xFITSImage(file_path)
        plt.figure('Toy point source')
        image.plot()
        ra0, dec0 = image.center()
        dra, ddec = image.delta()
        # Random sampling with no randomization---we expect (45., 45.)
        ra, dec = image.rvs_coordinates(1, randomize=False)
        self.assertAlmostEqual(ra, ra0)
        self.assertAlmostEqual(dec, dec0)
        # Random sampling with randomization---check bounds.
        ra, dec = image.rvs_coordinates(10000, randomize=True)
        self.assertAlmostEqual(max(ra - ra0), -0.5 * dra, places=6)
        self.assertAlmostEqual(min(ra - ra0), 0.5 * dra, places=6)
        self.assertAlmostEqual(max(dec - dec0), 0.5 * ddec, places=6)
        self.assertAlmostEqual(min(dec - dec0), -0.5 * ddec, places=6)

    @unittest.skip("Crab model removed for the public release")
    def test_crab(self):
        """Run a quick simulation of the Crab nebula and create a finely
        binned count map.
        """
        pipeline.reset('crab_nebula', overwrite=True)
        evt_file_list = pipeline.xpobssim(duration=1000)
        cmap_file_list = pipeline.xpbin(*evt_file_list, algorithm='CMAP',
                                        mc=True, pixsize=1.)
        if sys.flags.interactive:
            plt.figure('Crab binned image test')
            pipeline.xpbinview(*cmap_file_list)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

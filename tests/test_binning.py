#!/usr/bin/env python
#
# Copyright (C) 2018--2022, the ixpeobssim team.
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


"""Unit test for the binning module.
"""

import unittest
import os

import numpy
numpy.random.seed(666)

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.utils.misc import pairwise
from ixpeobssim.binning.base import xEventBinningBase
from ixpeobssim.binning.misc import xBinnedMap
from ixpeobssim.binning.polarization import xBinnedPolarizationMapCube
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.utils.logging_ import logger



class TestBinning(unittest.TestCase):

    """Unit test for the binning module.
    """

    def test_make_binning(self):
        """Basic test of the binning-creation functions.
        """
        mkbin = xEventBinningBase.make_binning
        binning = mkbin('LIN', min_val=1., max_val=10., num_bins=10)
        print(binning)
        binning = mkbin('LOG', min_val=1., max_val=10., num_bins=10)
        print(binning)
        data = numpy.random.uniform(1., 10., size=10000)
        binning = mkbin('EQP', num_bins=10, bin_data=data)
        print(binning)
        binning = mkbin('EQP', num_bins=10, bin_data=data, min_val=1.,
                         max_val=10.)
        print(binning)
        binning = mkbin('LIST', bin_list=[1, 2, 3])
        print(binning)

    def test_pairwise(self):
        """
        """
        binning = [1, 2, 3, 4, 5]
        for min_, max_ in pairwise(binning):
            print(min_, max_)

    def test_pairwise_vector(self):
        """
        """
        binning = numpy.linspace(1, 5, 5)
        for min_, max_ in pairwise(binning):
            print(min_, max_)

    def test_enumerate_pairwise(self):
        """
        """
        binning = [1, 2, 3, 4, 5]
        for i, (min_, max_) in enumerate(pairwise(binning)):
            print(i, min_, max_)

    def test_write(self):
        """Test writing the combined polarization map for the three DUs.
        """
        pipeline.reset('toy_disk', overwrite=True)
        file_list = pipeline.xpobssim(duration=1000)
        file_list = pipeline.xpbin(*file_list, algorithm='PMAP', npix=3)
        pmap = xBinnedPolarizationMapCube.from_file_list(file_list)
        file_path = os.path.join(IXPEOBSSIM_DATA, 'toy_disk_pmap_combined.fits')
        pmap.write(file_path)
        pmap_combined = xBinnedPolarizationMapCube(file_path)
        for ext_name in xBinnedPolarizationMapCube.EXTENSIONS_NAMES:
            logger.info('Comparing extension %s', ext_name)
            a1 = pmap.__getattr__(ext_name)[0]
            a2 = pmap_combined.__getattr__(ext_name)[0]
            delta = numpy.full(a1.shape, 0.)
            mask = a1 > 0.
            delta[mask] = abs((a1[mask] - a2[mask]) / a1[mask])
            self.assertTrue(numpy.allclose(a1, a2))
        file_list = pipeline.xpbin(*pipeline.file_list(), algorithm='CMAP', npix=3)
        cmap = xBinnedMap.from_file_list(file_list)
        file_path = os.path.join(IXPEOBSSIM_DATA, 'toy_disk_cmap_combined.fits')
        cmap.write(file_path)



if __name__ == '__main__':
    unittest.main()

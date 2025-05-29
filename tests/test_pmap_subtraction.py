#!/usr/bin/env python
#
# Copyright (C) 2025, the ixpeobssim team.
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

"""Unit tests for pmap rescaling and subtration.
"""

from __future__ import print_function, division

import unittest

import numpy
from astropy.io import fits

from ixpeobssim.binning.polarization import xBinnedPolarizationMapCube
from ixpeobssim.core import pipeline

# pylint: disable=invalid-name
# pylint: disable=no-member


class TestPmapSubtraction(unittest.TestCase):

    """Unit test for pmap subtraction and rescaling.
    """
    @classmethod
    def simulate(cls):
        """Run the simulation
        """
        DURATION_SRC=50000
        DURATION_SRC_ONLY = 25000
        DURATION_BKG=100000
        pipeline.reset('toy_point_source_bkg', overwrite=True)
        cls.source_list = pipeline.xpobssim(duration=DURATION_SRC,
                                             saa=False, occult=False)
        pipeline.reset('instrumental_bkg_smcx1', overwrite=True)
        cls.bkg_list = pipeline.xpobssim(duration=DURATION_BKG,
                                          saa=False, occult=False)
        pipeline.reset('toy_point_source_subtest', overwrite=True)
        cls.src_only_list = pipeline.xpobssim(duration=DURATION_SRC_ONLY,
                                               saa=False, occult=False)
        cls.src_lt = numpy.array([fits.open(file)[0].header['LIVETIME'] for
                          file in cls.source_list]).sum()
        cls.bkg_lt = numpy.array([fits.open(file)[0].header['LIVETIME'] for
                          file in cls.bkg_list]).sum()
        cls.src_only_lt = numpy.array([fits.open(file)[0].header['LIVETIME'] for
                               file in cls.src_only_list]).sum()

    def test(self):
        """
        """
        self.simulate()
        bkg_rescale = self.src_lt / self.bkg_lt
        src_only_rescale = self.src_lt / self.src_only_lt
        #Create pmaps
        src_map = xBinnedPolarizationMapCube.from_file_list(pipeline.xpbin
                                            (alg='PMAP', *self.source_list))
        bkg_map = xBinnedPolarizationMapCube.from_file_list(pipeline.xpbin
                                            (alg='PMAP', *self.bkg_list))
        src_only_map = xBinnedPolarizationMapCube.from_file_list(pipeline.xpbin
                                            (alg='PMAP', *self.src_only_list))
        #Rescale everything to the source + background map
        bkg_map *= bkg_rescale
        src_only_map *= src_only_rescale
        #Apply subtractions
        src_map -= bkg_map
        src_map -= src_only_map
        #Check that the results are compatible with zero (3 sigma threshold)
        npix = src_map.map_shape()[0]*src_map.map_shape()[1]
        assert(numpy.abs(src_map.hdu_list['I'].data.sum())<
               npix*3*src_map.hdu_list['I_ERR'].data.sum())

if __name__ == '__main__':
    unittest.main()

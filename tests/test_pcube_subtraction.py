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

"""Unit tests for pcube rescaling and subtration.
"""

from __future__ import print_function, division

import unittest

import numpy
from astropy.io import fits

from ixpeobssim.binning.polarization import xBinnedPolarizationCube
from ixpeobssim.core import pipeline
from ixpeobssim.instrument.gpd import GPD_DEFAULT_FIDUCIAL_HALF_SIDE_X,\
                                        GPD_DEFAULT_FIDUCIAL_HALF_SIDE_Y
from ixpeobssim.instrument.mma import fiducial_backscal

# pylint: disable=invalid-name
# pylint: disable=no-member

class TestPCUBESubtraction(unittest.TestCase):

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
        #Create source + background cubes
        src_reg = pipeline.xpselect(*self.source_list, rad=2.)
        src_cube_path = pipeline.xpbin(*src_reg, alg='PCUBE', emin=2.,
                                       emax=8., ebins=1)
        src_cube = xBinnedPolarizationCube.from_file_list(src_cube_path)
        #Create background cube
        bkg_cube = xBinnedPolarizationCube.from_file_list(pipeline.xpbin
                                    (alg='PCUBE', *self.bkg_list, ebins=1))
        #Create source only cube
        src_only_reg = pipeline.xpselect(*self.src_only_list, rad=2.)
        src_only_cubes_path = pipeline.xpbin(*src_only_reg, alg='PCUBE', emin=2.,
                                             emax=8., ebins=1)
        src_only_cube = xBinnedPolarizationCube.from_file_list(src_only_cubes_path)
        #Rescale everything to the source + background cube: first the livetime
        bkg_cube *= bkg_rescale
        src_only_cube *= src_only_rescale
        #Now also the backscal for the background
        full_backscal = fiducial_backscal(GPD_DEFAULT_FIDUCIAL_HALF_SIDE_X,
                                        GPD_DEFAULT_FIDUCIAL_HALF_SIDE_Y)
        bkg_cube *= src_cube.backscal()/full_backscal
        #Apply subtractions
        src_cube -= bkg_cube
        #Check that the results are compatible with zero (3 sigma threshold)
        assert(numpy.abs(src_cube.hdu_list[1].data['PD'][0]-
                         src_only_cube.hdu_list[1].data['PD'][0]) <=
                         3*numpy.abs(src_cube.hdu_list[1].data['PD_ERR'] +
                         src_only_cube.hdu_list[1].data['PD_ERR']))

if __name__ == '__main__':
    unittest.main()

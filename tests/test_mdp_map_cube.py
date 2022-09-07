#!/usr/bin/env python
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

import numpy

import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.binning.polarization import xBinnedMDPMapCube
from ixpeobssim.evt.kislat2015 import xStokesAnalysis
from ixpeobssim.utils.logging_ import logger



class TestMDPMapCube(unittest.TestCase):

    """Unit test for mdp map.

    """

    @classmethod
    def setUpClass(cls):
        """
        """
        pipeline.reset('toy_disk', overwrite=True)
        file_list = pipeline.xpobssim(duration=10000, seed=1)
        cls.file_list = pipeline.xpbin(*file_list, algorithm='MDPMAPCUBE',
            ebins=1, emin=2., emax=8., npix=5)

    @unittest.skip('Need to understand what we are doing with bins with no counts')
    def test_nan(self):
        """
        """
        for file_path in self.file_list:
            mdp_cube = xBinnedMDPMapCube(file_path)
            print(mdp_cube.I, mdp_cube.MU)
        mdp_cube = xBinnedMDPMapCube.from_file_list(file_list)
        print(mdp_cube.COUNTS, mdp_cube.MU, mdp_cube.I, mdp_cube.W2)

    def test_mdp_map_cube(self):
        """Simple method to simulate a toy disk source and bin an
        mdp_cube and test the values of the mdp.
        """
        mdp_cube = xBinnedMDPMapCube.from_file_list(self.file_list)
        mask = numpy.isfinite(mdp_cube.MU)
        mu = numpy.mean(mdp_cube.MU[mask])
        total_counts = numpy.sum(mdp_cube.COUNTS[mask])
        I = numpy.sum(mdp_cube.I[mask])
        W2 = numpy.sum(mdp_cube.W2[mask])
        overall_mdp = xStokesAnalysis.calculate_mdp99(mu, I, W2)
        logger.info('Total counts = %d, overall MDP = %.3f', total_counts, overall_mdp)
        pixel = 0, 2, 2
        mdp = mdp_cube.MDP_99[pixel]
        cnts_bin = mdp_cube.COUNTS[pixel]
        mu_bin = mdp_cube.MU[pixel]
        scale_mu = mu_bin / mu
        scale_factor = numpy.sqrt(cnts_bin / total_counts) * scale_mu
        scaled_mdp = mdp * scale_factor
        delta = (overall_mdp - scaled_mdp) / scaled_mdp
        self.assertTrue(abs(delta) < 0.1, 'delta = %f' % delta)



if __name__ == '__main__':
    unittest.main()

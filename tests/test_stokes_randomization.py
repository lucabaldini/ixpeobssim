#!/usr/bin/env python
#
# Copyright (C) 2022, the ixpeobssim team.
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

from ixpeobssim.binning.polarization import xBinnedPolarizationCube
import ixpeobssim.config.toy_point_source as toy_point_source
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.logging_ import logger


class TestStokesRandomization(unittest.TestCase):

    """Unit test for the randomization of the Stokes parameters.
    """

    @classmethod
    def setUpClass(cls):
        """
        """
        pipeline.reset('toy_point_source', overwrite=True)
        file_list = pipeline.xpobssim(duration=1000.)
        cls.evt_file_list = pipeline.xpselect(*file_list, emin=2., emax=8.)
        cls.pcube_file_list = pipeline.xpbin(*cls.evt_file_list, algorithm='PCUBE', ebins=1)

    def test_randomizing(self):
        """Test the randomization of the Stokes parameters.

        We want to make sure that the randomized Stokes parameters are properly
        normalized, and that the full output is unpolarized.
        """
        file_list = pipeline.xpstokesrandom(*self.evt_file_list)
        # Open all the randomized files and verify that the Stokes parameters are
        # properly normalized.
        for file_path in file_list:
            event_file = xEventFile(file_path)
            q = event_file.q_data()
            u = event_file.u_data()
            self.assertTrue(numpy.allclose(q**2. + u**2., 4.))
        file_list = pipeline.xpbin(*file_list, algorithm='PCUBE', ebins=1)
        pcube_original = xBinnedPolarizationCube.from_file_list(self.pcube_file_list)
        pcube_random = xBinnedPolarizationCube.from_file_list(file_list)
        logger.info(pcube_original)
        logger.info(pcube_random)
        for field in ('E_MEAN', 'COUNTS', 'MU', 'W2', 'N_EFF', 'FRAC_W', 'MDP_99'):
            self.assertTrue(numpy.allclose(pcube_original.__getattr__(field),
                pcube_random.__getattr__(field)))
        self.assertTrue((pcube_random.SIGNIF < 5.).all())



if __name__ == '__main__':
    unittest.main()

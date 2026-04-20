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

from ixpeobssim.binning.fmt import xBinTableHDUPCUBE
from ixpeobssim.binning.polarization import xBinnedPolarizationCube, xBinnedPolarizationMapCube
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca

if sys.flags.interactive:
    plt.ion()


class TestPolarizationMapCube(unittest.TestCase):

    """Unit test for the binned polarization cubes.
    """

    @classmethod
    def setUpClass(cls):
        """Run a small simulation and create the polarization cubes.
        """
        pipeline.reset('toy_gauss_disk', overwrite=True)
        cls.file_list = pipeline.xpobssim(duration=10000., seed=13)

    def test_sum_pixels(self):
        """Test the polarization cube addition.
        """
        kwargs = dict(algorithm='PCUBE', ebins=1, mc=True)
        file_list = pipeline.xpbin(*self.file_list, **kwargs)
        pcube = xBinnedPolarizationCube.from_file_list(file_list)
        kwargs = dict(algorithm='PMAPCUBE', ebins=1, mc=True, npix=200)
        file_list = pipeline.xpbin(*self.file_list, **kwargs)
        pmap = xBinnedPolarizationMapCube.from_file_list(file_list)
        mask = pcube.I[0] > 0.

        for key in ['EFLUX','EFLUXERR']:
            val1 = pcube.__getattr__(key)
            val2 = pmap.__getattr__(key)
            if key in ['EFLUX']:
                logger.info('%s -> %s %s', key, val1, pmap._sum_image_pixels(val2, mask))
                self.assertTrue(numpy.allclose(val1, pmap._sum_image_pixels(val2, mask), rtol=1.e-3, atol=1e-11))
            if key in ['EFLUXERR']:
                logger.info('%s -> %s %s', key, val1, numpy.sqrt(numpy.sum(val2**2)) )
                self.assertTrue(numpy.allclose(val1, numpy.sqrt(numpy.sum(val2**2)), rtol=1.e-3, atol=1e-12))



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

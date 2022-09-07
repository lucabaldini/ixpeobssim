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

import ixpeobssim.config.toy_point_source as srcmod

if sys.flags.interactive:
    plt.ion()


class TestPolarizationMapCube(unittest.TestCase):

    """Unit test for the binned polarization cubes.
    """

    @classmethod
    def setUpClass(cls):
        """Run a small simulation and create the polarization cubes.
        """
        pipeline.reset('toy_point_source', overwrite=True)
        cls.file_list = pipeline.xpobssim(duration=10000., seed=13)

    def test_simplest(self):
        """Compare a polarization cube and a map cube with one bin, and make
        sure the values are identical.
        """
        kwargs = dict(algorithm='PCUBE', ebins=1, suffix='pcube_simplest')
        file_list = pipeline.xpbin(*self.file_list, **kwargs)
        pcube = xBinnedPolarizationCube.from_file_list(file_list)
        kwargs = dict(algorithm='PMAP', suffix='pmap_simplest', npix=1)
        file_list = pipeline.xpbin(*self.file_list, **kwargs)
        pmap = xBinnedPolarizationMapCube.from_file_list(file_list)
        for key in xBinTableHDUPCUBE.POL_COL_NAMES:
            val1 = pcube.__getattr__(key)
            val2 = pmap.__getattr__(key)
            logger.info('%s -> %s %s', key, val1, val2)
            self.assertTrue(numpy.allclose(val1, val2[0, 0, 0]),
                '%s, delta = %s' % (key, val1 - val2[0, 0, 0]))
        pmap.plot()

    def test_energy(self, ebins=3):
        """Test with one spatial bin and multiple energy bins.
        """
        kwargs = dict(algorithm='PCUBE', ebins=ebins, suffix='pcube_energy')
        file_list = pipeline.xpbin(*self.file_list, **kwargs)
        pcube = xBinnedPolarizationCube.from_file_list(file_list)
        kwargs = dict(algorithm='PMAPCUBE', ebins=ebins, suffix='pmapcube_energy', npix=1)
        file_list = pipeline.xpbin(*self.file_list, **kwargs)
        pmap_cube = xBinnedPolarizationMapCube.from_file_list(file_list)
        for key in xBinTableHDUPCUBE.POL_COL_NAMES:
            val1 = pcube.__getattr__(key)
            val2 = pmap_cube.__getattr__(key)
            logger.info('%s -> %s %s', key, val1, val2[:, 0, 0])
            self.assertTrue(numpy.allclose(val1, val2[:, 0, 0]),
                '%s, delta = %s' % (key, val1 - val2[:, 0, 0]))

    def test_sky(self, npix=10):
        """
        """
        kwargs = dict(algorithm='PMAPCUBE', ebins=1, suffix='pmapcube_sky', npix=npix)
        file_list = pipeline.xpbin(*self.file_list, **kwargs)
        pmap_cube = xBinnedPolarizationMapCube.from_file_list(file_list)
        pmap_cube.plot()



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

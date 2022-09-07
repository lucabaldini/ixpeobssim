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

from ixpeobssim.binning.polarization import xBinnedPolarizationCube
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca

import ixpeobssim.config.toy_pollin as srcmod

if sys.flags.interactive:
    plt.ion()


class TestPolarizationCube(unittest.TestCase):

    """Unit test for the binned polarization cubes.
    """

    @classmethod
    def setUpClass(cls):
        """Run a small simulation and create the polarization cubes.
        """
        pipeline.reset('toy_pollin', overwrite=True)
        file_list = pipeline.xpobssim(duration=10000., seed=13)
        cls.file_list = pipeline.xpbin(*file_list, algorithm='PCUBE', ebins=3)

    def test_plot(self):
        """Plot the thing.
        """
        cube = xBinnedPolarizationCube.from_file_list(self.file_list)
        energy = numpy.linspace(2., 8., 20)
        plt.figure('Polarization degree')
        cube.plot_polarization_degree()
        plt.plot(energy, srcmod.pol_deg(energy))
        setup_gca(grids=True)
        plt.figure('Polarization angle')
        cube.plot_polarization_angle()
        plt.plot(energy, numpy.degrees(srcmod.pol_ang(energy)))
        setup_gca(grids=True)

    def test_addition(self):
        """Test the polarization cube addition.
        """
        file_path = self.file_list[0]
        cube1 = xBinnedPolarizationCube(file_path)
        cube2 = xBinnedPolarizationCube(file_path)
        cube2 += cube1
        for key in ['ENERG_LO', 'ENERG_HI', 'E_MEAN', 'MU', 'COUNTS', 'W2', 'I',
            'Q', 'U', 'PD', 'PD_ERR', 'PA', 'PA_ERR']:
            val1 = cube1.__getattr__(key)
            val2 = cube2.__getattr__(key)
            logger.info('%s\n  -> %s\n     %s', key, val1, val2)
            if key in ['COUNTS', 'W2', 'I', 'Q', 'U']:
                self.assertTrue(numpy.allclose(2 * val1, val2))
            elif key in ['PD_ERR', 'PA_ERR']:
                self.assertTrue(numpy.allclose(val1, numpy.sqrt(2.) * val2, rtol=1.e-3))
            else:
                self.assertTrue(numpy.allclose(val1, val2))



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

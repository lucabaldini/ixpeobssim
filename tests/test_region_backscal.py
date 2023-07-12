#!/usr/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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

"""Unit test for the BACKSCAL keyword support in region-selected sources.
"""

from __future__ import print_function, division

import unittest
import sys
import os

import numpy
from ixpeobssim.core import pipeline
from ixpeobssim.evt.event import xEventFile
from ixpeobssim import IXPEOBSSIM_CONFIG, IXPEOBSSIM_DATA, IXPEOBSSIM_TEST_DATA
from ixpeobssim.utils.logging_ import logger


class TestRegionBackscal(unittest.TestCase):

    """Unit test for the BACKSCAL keyword for xpselect with reg files
    """

    @classmethod
    def setUpClass(cls):
        """Test setup.
        """
        config_file_path = os.path.join(IXPEOBSSIM_CONFIG, 'toy_point_source.py')
        sim_file_path = os.path.join(IXPEOBSSIM_DATA, 'backscal_test')
        logger.info("Generating a toy point source from %s for selection testing...",
            config_file_path)
        cls.file_list = pipeline.xpobssim(configfile = config_file_path, duration=1000,
            outfile=sim_file_path)

    def test_region_backscal(self):
        """Compare the reg file BACKSCAL with that of the ordinary circle cut of 1 arcmin.
        """
        evt_file_path = self.file_list[0]
        region_file_path = os.path.join(IXPEOBSSIM_TEST_DATA, 'test_reg_backscal.reg')
        logger.info('Selecting a circular region from %s...', region_file_path)
        reg_selected = pipeline.xpselect(evt_file_path, regfile=region_file_path, suffix='reg',
            overwrite=True)
        reg_backscal = xEventFile(*reg_selected).backscal()
        logger.info('Selecting a circular region with built-in xpselect function...')
        circ_selected = pipeline.xpselect(evt_file_path, rad=1, suffix='circ', overwrite=True)
        circ_backscal = xEventFile(*circ_selected).backscal()
        logger.info('Region backscal: %.3f arcsec^2', reg_backscal)
        logger.info('Circle backscal: %.3f arcsec^2', circ_backscal)
        logger.info('Ratio: %.6f', reg_backscal / circ_backscal)
        assert numpy.allclose(reg_backscal, circ_backscal, rtol=1e-2)

    def test_region_backscal_inverse(self):
        """Compare the reg file BACKSCAL (with the --reginv flag set) with a
        simple annulus with a 1 arcmin inner radius.
        """
        evt_file_path = self.file_list[0]
        region_file_path = os.path.join(IXPEOBSSIM_TEST_DATA, 'test_reg_backscal.reg')
        logger.info('Selecting an inverted circular region from %s...', region_file_path)
        reg_selected = pipeline.xpselect(evt_file_path, regfile=region_file_path,
            suffix='reginv', reginv=True, overwrite=True)
        reg_backscal = xEventFile(*reg_selected).backscal()
        logger.info('Selecting an annulus with built-in xpselect function...')
        circ_selected = pipeline.xpselect(evt_file_path, innerrad=1, suffix='annulus',
            overwrite=True)
        circ_backscal = xEventFile(*circ_selected).backscal()
        logger.info('Region backscal: %.3f arcsec^2', reg_backscal)
        logger.info('Circle backscal: %.3f arcsec^2', circ_backscal)
        logger.info('Ratio: %.6f', reg_backscal / circ_backscal)
        assert numpy.allclose(reg_backscal, circ_backscal, rtol=1e-2)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

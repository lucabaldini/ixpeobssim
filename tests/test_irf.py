#!/usr/bin/env python
#
# Copyright (C) 2018, the ixpeobssim team.
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
import os

from ixpeobssim.instrument import DU_IDS
from ixpeobssim.irf import DEFAULT_IRF_NAME, load_irf_set, load_arf, load_mrf
from ixpeobssim.irf.caldb  import IRF_TYPES
from ixpeobssim.irf.caldb  import irf_file_name, irf_folder_path, irf_file_path, parse_irf_name
from ixpeobssim.utils.logging_ import logger


class TestIrf(unittest.TestCase):

    """Unit test for the irf __init__.py module.
    """

    def test_folder_path(self):
        """Make sure we get the path to the irf folder right.
        """
        for irf_type in IRF_TYPES:
            # The GPD quantum efficiency files are not shipped with ixpeobssim
            # and therefore we skip the corresponding tests, here.
            if irf_type in ('qe', ):
                continue
            folder_path = irf_folder_path(irf_type)
            logger.info('%s -> %s', irf_type, folder_path)
            self.assertTrue(os.path.isdir(folder_path))

    def test_file_name(self, base='ixpe_', du_id=1, version=1):
        """
        For reference, here is the naming conventions for the actual CALDB

        ixpe_d1_20170101_01.arf
        ixpe_d1_20170101_alpha075_01.arf
        ixpe_d1_20170101_01.mrf
        ixpe_d1_20170101_alpha075_01.mrf
        ixpe_d1_20170101_01.rmf
        ixpe_d1_20170101_alpha075_01.rmf
        ixpe_d1_20170101_mfact_01.fits
        ixpe_d1_20170101_mfact_alpha075_01.fits
        """
        for irf_type in IRF_TYPES:
            for intent in ('_obssim_', '_obssim_alpha075_'):
                file_name = irf_file_name(base, du_id, irf_type, intent, version)
                print(file_name)

    def test_version_switch(self, base='ixpe', intent='obssim', du_id=1):
        """Test the version switch between versions 9 and 10.
        """
        for irf_type in IRF_TYPES:
            for version in (9, 10):
                file_name = irf_file_name(base, du_id, irf_type, intent, version)
                print(file_name)

    def test_irf_name(self):
        """Test the mechanisms for splitting IRF names.
        """
        for irf_name in ('ixpe:legacy_stdcut:v6', 'ixpe:legacy_stdcut:v6'):
            logger.info('%s -> %s', irf_name, parse_irf_name(irf_name))

    def test_irf_file_path(self, irf_name=DEFAULT_IRF_NAME):
        """Make sure that the actual IRF file exist.
        """
        logger.info('Checking files on disk for IRF name %s', irf_name)
        for irf_type in IRF_TYPES:
            # The GPD quantum efficiency files are not shipped with ixpeobssim
            # and therefore we skip the corresponding tests, here.
            if irf_type in ('qe', ):
                continue
            for du_id in [1, 2, 3]:
                file_path = irf_file_path(irf_name, du_id, irf_type, check_file=True)
                logger.info(file_path)
                self.assertTrue(os.path.isfile(file_path))

    def test_irf_load(self, irf_name=DEFAULT_IRF_NAME):
        """
        """
        for du_id in DU_IDS:
            irf_set = load_irf_set(irf_name, du_id)
            for irf in (irf_set.aeff, irf_set.vign, irf_set.edisp, irf_set.psf, irf_set.modf, irf_set.mrf):
                self.assertEqual(irf.du_id, du_id)
                logger.info(irf.header_comments())
                logger.info(irf.file_path)
                logger.info(irf.irf_type)

    def test_simple_weights(self):
        """Quick test for loading the arf and mrf files with the SIMPLE weighting
        prescription.
        """
        for du_id in DU_IDS:
            # Load the actual weighted response files with the SIMPLE weighting scheme.
            aeff = load_arf('ixpe:obssim_alpha075:v12', du_id, simple_weighting=True)
            self.assertTrue('simple' in aeff.file_path)
            mrf = load_mrf('ixpe:obssim_alpha075:v12', du_id, simple_weighting=True)
            self.assertTrue('simple' in mrf.file_path)
            # Make sure that asking for the simple weighting prescription in the
            # unweighted case is rising a RuntimeError.
            self.assertRaises(RuntimeError, load_arf, 'ixpe:obssim:v12', du_id, simple_weighting=True)
            self.assertRaises(RuntimeError, load_mrf, 'ixpe:obssim:v12', du_id, simple_weighting=True)
            # Load a full IRF set.
            irf_set = load_irf_set('ixpe:obssim_alpha075:v12', du_id, simple_weighting=True)
            self.assertTrue('simple' in irf_set.aeff.file_path)
            self.assertTrue('simple' in irf_set.mrf.file_path)

    def test_gray_filter(self):
        """Quick test for the response files with the gray filter.
        """
        for du_id in DU_IDS:
            irf_name = 'ixpe:obssim:v12'
            aeff = load_arf(irf_name, du_id, gray_filter=True)
            self.assertTrue('gray' in aeff.file_path)
            mrf = load_arf(irf_name, du_id, gray_filter=True)
            self.assertTrue('gray' in mrf.file_path)
            irf_name = 'ixpe:obssim_alpha075:v12'
            aeff = load_arf(irf_name, du_id, gray_filter=True)
            self.assertTrue('gray' in aeff.file_path)
            mrf = load_arf(irf_name, du_id, gray_filter=True)
            self.assertTrue('gray' in mrf.file_path)
            aeff = load_arf(irf_name, du_id, simple_weighting=True, gray_filter=True)
            self.assertTrue('gray' in aeff.file_path)
            self.assertTrue('simple' in aeff.file_path)
            mrf = load_arf(irf_name, du_id, simple_weighting=True, gray_filter=True)
            self.assertTrue('gray' in mrf.file_path)
            self.assertTrue('simple' in mrf.file_path)



if __name__ == '__main__':
    unittest.main()

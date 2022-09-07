#!/usr/bin/env python
#
# Copyright (C) 2019, the ixpeobssim team.
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


"""Unit test for reading ds9 region files.
"""


import unittest
import os

from ixpeobssim import IXPEOBSSIM_CONFIG_REG
from ixpeobssim.utils.astro import read_ds9


class TestDs9(unittest.TestCase):

    """Unit test for ds9 regions.
    """

    def test_casa(self):
        """
        """
        file_path = os.path.join(IXPEOBSSIM_CONFIG_REG, 'casa_multiple.reg')
        region_list = read_ds9(file_path)
        for region in region_list:
            print(region)

    def test_cena(self):
        """
        """
        file_path = os.path.join(IXPEOBSSIM_CONFIG_REG, 'cena_jet+core.reg')
        region_list = read_ds9(file_path)
        for region in region_list:
            print(region)



if __name__ == '__main__':
    unittest.main()

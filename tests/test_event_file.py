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

import numpy
import unittest

import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.evt.event import xEventFile



class TestEventFile(unittest.TestCase):

    """.
    """

    @classmethod
    def setUpClass(cls):
        """
        """
        pipeline.reset('toy_casa', overwrite=True)
        cls.file_list = pipeline.xpobssim(duration=1000.)

    def test_wcs(self):
        """
        """
        for file_path in self.file_list:
            evt_file = xEventFile(file_path)
            ra0, dec0 = evt_file.wcs_reference()
            self.assertAlmostEqual(ra0, evt_file.primary_header.get('RA_OBJ'))
            self.assertAlmostEqual(ra0, evt_file.primary_header.get('RA_PNT'))
            self.assertAlmostEqual(dec0, evt_file.primary_header.get('DEC_OBJ'))
            self.assertAlmostEqual(dec0, evt_file.primary_header.get('DEC_PNT'))



if __name__ == '__main__':
    unittest.main()

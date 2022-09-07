#!/usr/bin/env python
#
# Copyright (C) 2016--2018, the ixpeobssim team.
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

import numpy

import ixpeobssim.core.pipeline as pipeline
from ixpeobssim import IXPEOBSSIM_CONFIG_REG
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.astro import angular_separation, read_ds9
from ixpeobssim.utils.units_ import degrees_to_arcmin
from ixpeobssim.utils.matplotlib_ import plt


class TestSelect(unittest.TestCase):

    """Unit test for simulation and analysis pipeline(s).
    """

    def test_toy_point_source(self):
        """Simulate a simple source with a power-law spectrum.
        """
        pipeline.reset('toy_point_source', overwrite=True)
        file_list = pipeline.xpobssim(duration=1000., saa=False, occult=False)
        pipeline.xpselect(*file_list, emin=5, emax=6)
        pipeline.xpselect(*file_list, rad=0.1)

    def test_circular_selection(self):
        """This was added to investigate issue #282.

        Basically we are generating an extended disk, cutting out the inner
        1-arcmin circle based on a region file, and seeing what happens.
        """
        reg_file_path = os.path.join(IXPEOBSSIM_CONFIG_REG, 'toy_disc_circle.reg')
        # Retrieve the region information from the region file.
        circle = read_ds9(reg_file_path)[0]
        ra0, dec0, rad = circle.center.ra.deg, circle.center.dec.deg, circle.radius
        rad = degrees_to_arcmin(rad)
        # Run the simulation.
        pipeline.reset('toy_disk', overwrite=True)
        file_list = pipeline.xpobssim(duration=100000.)
        file_list = pipeline.xpselect(*file_list, regfile=reg_file_path)
        for file_path in file_list:
            event_file = xEventFile(file_path)
            ra, dec = event_file.sky_position_data()
            cdelt1, cdelt2 = event_file._wcs.wcs.cdelt
            delta = degrees_to_arcmin(max(abs(cdelt1), abs(cdelt2)))
            dist = angular_separation(ra, dec, ra0, dec0)
        pipeline.xpbin(*file_list, algorith='CMAP')



if __name__ == '__main__':
    unittest.main()

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


import ixpeobssim.core.pipeline as pipeline
from ixpeobssim import IXPEOBSSIM_CONFIG
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.profile import xMemoryProfiler
from ixpeobssim.evt.event import xEventFile

profiler = xMemoryProfiler()

class TestMemoryProfile(unittest.TestCase):

    """
    """

    def test_xpselect(self):
        """
        """
        pipeline.reset('toy_point_source')
        logger.info(profiler)
        file_list = pipeline.xpobssim(duration=10000.)
        logger.info(profiler)
        samples = []
        for i in range(5):
            for file_path in file_list:
                pipeline.xpselect(file_path, suffix='sel%d' % i, overwrite=True,
                                  emin=4.)
                logger.info(profiler)
                samples.append(profiler.last_avail_mem)
        logger.info(profiler)
        logger.info(samples)


    def test_open_event_file(self):
        """
        """
        pipeline.reset('toy_point_source')
        logger.info(profiler)
        file_list = pipeline.xpobssim(duration=10000.)
        for i in range(5):
            for file_path in file_list:
                event_file = xEventFile(file_path)
                ra, dec = event_file.sky_position_data()
                logger.info(profiler)



if __name__ == '__main__':
    unittest.main()
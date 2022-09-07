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
import numpy

import ixpeobssim.core.pipeline as pipeline

from ixpeobssim.evt.event import xEventFile


class TestTriggerId(unittest.TestCase):

    """
    """

    def base_test(self, source_name):
        """
        """
        pipeline.reset(source_name, overwrite=True)
        file_list = pipeline.xpobssim(duration=100., saa=False, occult=False)
        for file_path in file_list:
            event_file = xEventFile(file_path)
            trg_id = event_file.trigger_id_data()
            target = numpy.arange(1, event_file.num_events() + 1)
            self.assertTrue(numpy.array_equal(trg_id, target))

    def test_toy_point_source(self):
        """
        """
        self.base_test('toy_point_source')

    def test_multiple_sources(self):
        """
        """
        self.base_test('toy_multiple_sources')

    def test_toy_periodic_source(self):
        """
        """
        self.base_test('toy_periodic_source')



if __name__ == '__main__':
    unittest.main()

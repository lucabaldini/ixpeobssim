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

from __future__ import print_function, division

import os
import sys
import unittest

import numpy
import matplotlib.pyplot as plt

from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.srcmodel.ephemeris import xOrbitalEphemeris
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.time_ import MISSION_START_UNIX_TIME
import ixpeobssim.core.pipeline as pipeline

if sys.flags.interactive:
    plt.ion()



class TestXpphase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Test initialization.
        """
        cls.file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'parfiles', 'SAXJ1808.4-3658.par')
        cls.ephem = xOrbitalEphemeris.from_file(cls.file_path)
        cls.obs_length = 1000
        logger.info(cls.ephem)

    def test_get_phase(self):
        """Basic test
        """
        start_met = MISSION_START_UNIX_TIME
        t = numpy.linspace(0, self.obs_length, 100)
        ph_array = self.ephem.fold(t, start_met)
        plt.figure()
        plt.plot(t, ph_array)

    def test_xpphase(self):
        """
        """
        pipeline.reset('toy_disk', overwrite=True)
        file_list = pipeline.xpobssim(duration=10000)
        file_list = pipeline.xpphase(*file_list, parfile=self.file_path, suffix='folded')



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

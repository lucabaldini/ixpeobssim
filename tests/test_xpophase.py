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
from ixpeobssim.bin.xpophase import barycorr, get_ophase
from ixpeobssim.utils.logging_ import logger
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

    def test_barycorr(self):
        """Basic test
        """
        t = numpy.linspace(0, self.obs_length, 100)
        t = barycorr(t, self.ephem, ra=45, dec=45)

    def test_get_ophase(self):
        """Basic test
        """
        t = numpy.linspace(0, self.obs_length, 100)
        ph_array = get_ophase(t, self.ephem)
        plt.figure()
        plt.plot(t, ph_array)

    def test_xpophase(self):
        """
        """
        pipeline.reset('toy_disk', overwrite=True)
        file_list = pipeline.xpobssim(duration=10000)
        file_list = pipeline.xpophase(*file_list, parfile=self.file_path,
                                      bary=True, ra=45, dec=45, suffix='folded')


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

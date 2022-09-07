#!/usr/bin/env python
#
# Copyright (C) 2018--2022, the ixpeobssim team.
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

from ixpeobssim import IXPEOBSSIM_BIN
import sys
sys.path.append(IXPEOBSSIM_BIN)
sys.dont_write_bytecode = 1
from xpbin import BIN_ALGS
from ixpeobssim.binning import read_binned_file_list, BINNING_READ_DICT
from ixpeobssim.utils.matplotlib_ import plt
import ixpeobssim.core.pipeline as pipeline

if sys.flags.interactive:
    plt.ion()



class TestBin(unittest.TestCase):

    """Unit test for simulation and analysis pipeline(s).
    """

    def test_toy_point_source(self):
        """Run a quick simulation and create a finely binned count map.
        """
        pipeline.reset('toy_point_source', overwrite=True)
        evt_file_list = pipeline.xpobssim(duration=1000)
        cmap_file_list = pipeline.xpbin(*evt_file_list, algorithm='CMAP', mc=True)
        if sys.flags.interactive:
            plt.figure('Toy point source binned image test')
            pipeline.xpbinview(*cmap_file_list)

    @unittest.skip("Disengaged for the public release")
    def test_crab_pulsar(self):
        """Simulate a simple periodic source.
        """
        pipeline.reset('crab_pulsar', overwrite=True)
        file_list = pipeline.xpobssim(duration=2000., saa=False, occult=False)
        file_list = pipeline.xpphase(*file_list, **crab_pulsar.ephemeris.dict())
        for alg in BIN_ALGS:
            bin_file_list = pipeline.xpbin(*file_list, algorithm=alg)
            if alg in BINNING_READ_DICT:
                _sum = read_binned_file_list(alg, bin_file_list)
                _sum.plot()




if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

#!/usr/bin/env python
#
# Copyright (C) 2021--2022, the ixpeobssim team.
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

import sys
import unittest

from ixpeobssim.bin.xpbin import BIN_ALGS
from ixpeobssim.binning.misc import xBinnedMap
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.irf import DEFAULT_IRF_NAME
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt

import ixpeobssim.config.toy_periodic_source as input_model

if sys.flags.interactive:
    plt.ion()
from ixpeobssim.utils.environment import PYXSPEC_INSTALLED
if PYXSPEC_INSTALLED:
    import ixpeobssim.evt.xspec_ as xspec_



class TestStripmc(unittest.TestCase):

    """
    """

    @classmethod
    def setUpClass(self):
        """Run a short simulation and strip the Monte Carlo information.
        """
        pipeline.reset('toy_periodic_source', overwrite=True)
        file_list = pipeline.xpobssim(duration=10000.)
        file_list = pipeline.xpstripmc(*file_list)
        self.file_list = pipeline.xpphase(*file_list, **input_model.ephemeris.dict())

    def test_xpbin(self):
        """Test xpbin in all the possible flavors.
        """
        for alg in BIN_ALGS:
            logger.info('Testing xpbin in the "--%s" flavor...', alg)
            file_list = pipeline.xpbin(*self.file_list, algorithm=alg, irfname=DEFAULT_IRF_NAME)

    def test_xpselect(self):
        """Test xpselect.
        """
        pipeline.xpselect(*self.file_list, emin=2.)
        pipeline.xpselect(*self.file_list, rad=2.)

    @unittest.skip('Needs to be done properly with a non-period source')
    def test_xpxspec(self):
        """Test xpxspec.
        """
        if not PYXSPEC_INSTALLED:
            return
        file_list = pipeline.file_list('nomc', 'pha1*')
        fit_output = pipeline.xpxspec(*file_list, model='powerlaw * polconst', plot=sys.flags.interactive)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

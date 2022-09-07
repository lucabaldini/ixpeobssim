#!/usr/bin/env python
#
# Copyright (C) 2016, the ixpeobssim team.
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


"""Unit test for IRF plotting.
"""

import unittest
import sys

from ixpeobssim.irf import load_irf_set, DEFAULT_IRF_NAME
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()


class TestIrfPlot(unittest.TestCase):

    """
    """

    @classmethod
    def setUpClass(cls):
        """Setup.
        """
        cls.irf_name = DEFAULT_IRF_NAME
        cls.irf_set = load_irf_set(cls.irf_name, du_id=1)

    def test_irfplot(self):
        """Plot all the instrument response functions.
        """
        self.irf_set.aeff.plot()
        self.irf_set.vign.plot()
        self.irf_set.psf.plot()
        self.irf_set.modf.plot()
        self.irf_set.edisp.plot()



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

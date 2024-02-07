#!/usr/bin/env python
#
# Copyright (C) 2016--2019, the ixpeobssim team.
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
import sys

import numpy

from ixpeobssim.irf import DEFAULT_IRF_NAME
from ixpeobssim.srcmodel.spectrum import xSmearingMatrix
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color

if sys.flags.interactive:
    plt.ion()



class TestSmearingMatrix(unittest.TestCase):

    """Unit test for the IXPE smearing matrix.
    """

    def test(self, irf_name=DEFAULT_IRF_NAME, du_id=1, spectral_index=2.):
        """
        """
        matrix_std = xSmearingMatrix(irf_name, du_id, spectral_index, gray_filter=False)
        matrix_gry = xSmearingMatrix(irf_name, du_id, spectral_index, gray_filter=True)
        for chan in (50, 75):
            plt.figure('Smearing matrix slice at channel %d' % chan)
            matrix_std.vslice(chan).plot(label='Standard')
            matrix_gry.vslice(chan).plot(label='Gray filter')
            plt.text(0.75, 0.75, 'Channel %d' % chan, transform=plt.gca().transAxes)
            setup_gca(xlabel='Energy [keV]', ylabel='pdf', legend=True,
                logy=True, grids=True, ymin=1e-8)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

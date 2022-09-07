#!/usr/bin/env python
#
# Copyright (C) 2015, the ixpeobssim team.
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


"""Unit test for the irf.psf module.
"""


import numpy
import sys
import os
import unittest

from ixpeobssim import IXPEOBSSIM_CONFIG_FITS
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.irf import load_psf, DEFAULT_IRF_NAME
from ixpeobssim.irf.psf import gauss_king
from ixpeobssim.core.fitsio import xFITSImageBase
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()


class TestIxpePsf(unittest.TestCase):

    """Unit test for the IXPE point-spread function.
    """

    def test_plot(self):
        """
        """
        file_name = 'tycho_4p1_6p1_keV.fits'
        file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, file_name)
        img = xFITSImageBase(file_path)
        plt.figure('PSF circle')
        fig = img.plot()
        psf = load_psf(DEFAULT_IRF_NAME)
        psf.draw_psf_circle(img, hpd_value=True)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

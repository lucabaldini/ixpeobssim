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


"""Unit test for the animation module.
"""


import unittest
import sys
import os

import numpy

from ixpeobssim import IXPEOBSSIM_CONFIG_FITS
from ixpeobssim.config.casa_animation import ROI
from ixpeobssim.evt.animate import xMovingCircle, xSkyAnimation
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.units_ import arcmin_to_degrees

if sys.flags.interactive:
    plt.ion()


class TestAnimate(unittest.TestCase):

    """
    """

    def test_moving_circle(self, tmin=0., tmax=20., num_events=2500):
        """
        """
        positions = (tmin, 0., 0.), (tmax, 5., 0.)
        roi = xMovingCircle(positions)
        t = numpy.random.uniform(tmin, tmax, size=num_events)
        t.sort()
        ra = numpy.random.uniform(-0.05, 0.05, size=num_events)
        dec = numpy.random.uniform(-0.05, 0.05, size=num_events)
        mask = roi.event_mask(t, ra, dec, 0., 0.)
        plt.figure('Test animation')
        plt.plot(ra, dec, 'o')
        plt.plot(ra[mask], dec[mask], 'o')
        plt.gca().set_aspect('equal')
        setup_gca(xlabel='R.A. [degrees]', ylabel='Dec. [degrees]', grids=True)

    def test_casa(self):
        """Test roi path with Cas A.
        """
        file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'casa_1p5_3p0_keV.fits')
        img = xSkyAnimation(file_path)
        delta_ra, delta_dec, _ = ROI(numpy.linspace(ROI.tmin, ROI.tmax, 100))
        ra = img.ra0 + arcmin_to_degrees(delta_ra)
        dec = img.dec0 +  arcmin_to_degrees(delta_dec)
        x, y = img.image.wcs.wcs_world2pix(ra, dec, 0)
        plt.figure('Cas A animation')
        img.plot()
        plt.plot(x, y, color='white')




if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

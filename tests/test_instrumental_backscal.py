# Copyright (C) 2022, the ixpeobssim team.
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

from ixpeobssim.core.hist import xGpdMap2d, xHistogram2d
from ixpeobssim.instrument.gpd import GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_Y,\
    within_fiducial_rectangle, GPD_DEFAULT_FIDUCIAL_HALF_SIDE_X, GPD_DEFAULT_FIDUCIAL_HALF_SIDE_Y
from ixpeobssim.instrument.mma import apply_dithering, gpd_to_sky, fiducial_backscal
from ixpeobssim.srcmodel.roi import xModelComponentBase
from ixpeobssim.utils.astro import angular_separation
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.units_ import degrees_to_arcmin, arcmin_to_arcsec


if sys.flags.interactive:
    plt.ion()



class TestInstrumentalBackscal(unittest.TestCase):

    """Unit test for the backscal calculation for the instrumenatal background.
    """

    def test(self, num_events=1000000, duration=100000., start_met=0.,
        dither_params=(1.6, 907., 101., 449.)):
        """
        """
        # Generate events uniformly over the detector surface.
        logger.info('Generating %d events on the active surface...', num_events)
        args = num_events, GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_Y
        detx, dety = xModelComponentBase.uniform_rectangle(*args)
        t = xModelComponentBase.uniform_time(num_events, start_met, duration)
        # Apply the fiducial cut.
        mask = within_fiducial_rectangle(detx, dety)
        fiducial_events = mask.sum()
        logger.info('Trimming to the fiducial area, %d events remaining...', fiducial_events)
        detx, dety, t = detx[mask], dety[mask], t[mask]
        # Apply dithering and project into sky coordinates.
        ra_pnt, dec_pnt = apply_dithering(t, 0., 0., dither_params)
        ra, dec = gpd_to_sky(detx, dety, t, ra_pnt, dec_pnt, du_id=1, roll_angle=0.)

        plt.figure('Events in instrument coordinates')
        hist = xGpdMap2d(100).fill(detx, dety)
        hist.plot()

        plt.figure('Events in sky coordinates')
        side = 0.15
        binning = numpy.linspace(-side, side, 100)
        hist = xHistogram2d(binning, binning).fill(ra, dec)
        hist.plot()
        setup_gca(xlabel='R. A. [degrees]', ylabel='Dec. [degrees]')

        full_backscal = fiducial_backscal(GPD_DEFAULT_FIDUCIAL_HALF_SIDE_X, GPD_DEFAULT_FIDUCIAL_HALF_SIDE_Y)
        radius = numpy.linspace(1., 10., 50)
        backscal_ratio = []
        backscal_ratio_err = []
        angsep = degrees_to_arcmin(angular_separation(ra, dec, 0., 0.))
        for r in radius:
            mask = angsep <= r
            n = mask.sum()
            frac_events = n / fiducial_events
            backscal = numpy.pi * arcmin_to_arcsec(r)**2.
            frac_backscal = backscal / full_backscal
            ratio = frac_events / frac_backscal
            ratio_err = frac_events * (1. - frac_events) / numpy.sqrt(n) / frac_backscal
            backscal_ratio.append(ratio)
            backscal_ratio_err.append(ratio_err)
        backscal_ratio = numpy.array(backscal_ratio)
        backscal_ratio_err = numpy.array(backscal_ratio_err)

        plt.figure('Backscal calculation')
        plt.errorbar(radius, backscal_ratio, backscal_ratio_err, fmt='o')
        setup_gca(ymin=0., ymax=1.1, grids=True, xlabel='Extraction radius [arcmin]',
            ylabel='Fraction of events / fractional backscal')



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

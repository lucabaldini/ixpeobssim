#!/usr/bin/env python
#
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

from ixpeobssim.instrument.mma import apply_dithering
from ixpeobssim.instrument.sc import pointing_splines
from ixpeobssim.instrument.traj import xObservationTimeline
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, residual_plot
from ixpeobssim.utils.time_ import LAUNCH_MET
from ixpeobssim.utils.units_ import degrees_to_arcsec

if sys.flags.interactive:
    plt.ion()



class TestTimeline(unittest.TestCase):

    """Unit test for the timeline code.
    """

    def test_sc_data(self, duration=1000., sc_interval=10., target_ra=180.,
        target_dec=30., saa=True, occult=True):
        """Create a timeline object and verify the SC_DATA granularity that is
        necessary to recover the pointing with the necessary degree of precision.
        """
        start_met = LAUNCH_MET
        stop_met = start_met + duration
        # Create the timeline and retrieve the spacecraft data.
        timeline = xObservationTimeline(start_met, stop_met, target_ra, target_dec, saa, occult)
        dither_params = (1.6, 907., 101., 449.)
        sc_data = timeline.sc_data(sc_interval, dither_params)
        # Calculate the proper dithering splines from the spacecraft data.
        ra_spline, dec_spline = pointing_splines(sc_data)
        # Recalculate the dithering pattern over a much finer grid.
        _met = numpy.linspace(start_met, stop_met, int(duration) * 100)
        _ra, _dec = apply_dithering(_met, target_ra, target_dec, dither_params)
        # Estimate the errors.
        logger.info('Estimating the spline approximation errors on the dithering...')
        delta_ra = degrees_to_arcsec(_ra - ra_spline(_met))
        delta_dec = degrees_to_arcsec(_dec - dec_spline(_met))
        max_ra_error = max(abs(delta_ra))
        max_dec_error = max(abs(delta_dec))
        logger.info('Maximum R.A. error: %.2f arcsec', max_ra_error)
        logger.info('Maximum Dec. error: %.2f arcsec', max_dec_error)

        ax1, ax2 = residual_plot('RA timeline')
        plt.plot(sc_data['MET'], sc_data['RA_PNT'], 'o', label='SC_DATA')
        plt.plot(_met, _ra, label='True')
        ra_spline.plot(label='Spline approx.')
        setup_gca(grids=True, ylabel='R.A. [deg]', legend=True)
        plt.sca(ax2)
        plt.plot(_met, delta_ra)
        setup_gca(grids=True, xlabel='MET [s]', ylabel='Delta R.A. [arcsec]')

        ax1, ax2 = residual_plot('Dec timeline')
        plt.plot(sc_data['MET'], sc_data['DEC_PNT'], 'o', label='SC_DATA')
        plt.plot(_met, _dec, label='True')
        dec_spline.plot(label='Spline approx.')
        setup_gca(grids=True, ylabel='Declination [deg]', legend=True)
        plt.sca(ax2)
        plt.plot(_met, delta_dec)
        setup_gca(grids=True, xlabel='MET [s]', ylabel='Delta Dec. [arcsec]')


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

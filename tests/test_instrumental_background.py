#!/usr/bin/env python
#
# Copyright (C) 2020, the ixpeobssim team.
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

import numpy
import astropy.io.fits as fits

from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.config.instrumental_bkg import bkg
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.instrument.gpd import GPD_PHYSICAL_AREA, fiducial_area
import ixpeobssim.core.pipeline as pipeline


if sys.flags.interactive:
    plt.ion()


class TestInstrumentalBackground(unittest.TestCase):

    """Unit test for the instrumental background.
    """

    def test(self, duration=1.e6):
        """Note this is run with zero deadtime, no SAA nor Earth occultation in
        order to be able to compute the normalization correctly.
        """
        pipeline.reset('instrumental_bkg', overwrite=True)
        evt_file_list = pipeline.xpobssim(duration=duration, deadtime=0., saa=False,
            occult=False, scdata=False)
        binning = numpy.linspace(1., 12., 100)
        bin_width = binning[1] - binning[0]
        hist = xHistogram1d(binning, xlabel='Energy [keV]')
        for file_path in evt_file_list:
            with fits.open(file_path) as hdu_list:
                data = hdu_list['MONTE_CARLO'].data
                hist.fill(data['MC_ENERGY'])
        scale = duration * bin_width * (fiducial_area() / 100.) * 3.
        x = hist.bin_centers(0)
        y = hist.content
        chisq = (((y - scale * bkg.photon_spectrum(x)) / numpy.sqrt(y))**2.).sum()
        ndof = len(y)
        delta = abs((chisq - ndof) / numpy.sqrt(2. * ndof))
        self.assertTrue(delta <= 5.)

        plt.figure('Background spectrum')
        hist *= 1. / scale
        hist.plot(label='Simulation output')
        energy = numpy.linspace(1.0, 15, 100)
        plt.plot(energy, bkg.photon_spectrum(energy), label='Model')
        setup_gca(logx=True, logy=True, ymin=3.e-6, ymax=5.e-3, grids=True, legend=True)

    def test_gtis(self, duration=1.e6):
        """And this is running the same thing with non trivial GTIs to avoid
        choking on bugs a la #482
        """
        pipeline.reset('instrumental_bkg', overwrite=True)
        evt_file_list = pipeline.xpobssim(duration=duration, saa=True, occult=True,
            scdata=False)




if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

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


"""Unit test for the 2d part of the irf.psf module.
"""

import os
import sys
import unittest

import numpy

from ixpeobssim import IXPEOBSSIM_IRFGEN_DATA
from ixpeobssim.core.hist import xHistogram1d, xHistogram2d
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.instrument.mma import FOCAL_LENGTH
from ixpeobssim.irf.caldb import irf_file_path
from ixpeobssim.irf.psf import xPointSpreadFunction2d
from ixpeobssim.irf import load_psf, DEFAULT_IRF_NAME
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.units_ import degrees_to_arcmin

if sys.flags.interactive:
    plt.ion()


class TestIxpePsf2d(unittest.TestCase):

    """Unit test for the IXPE point-spread function in the 2d version.
    """

    _IRF_NAME = 'ixpe:obssim:v12'

    @unittest.skip('Waiting for the PSF files.')
    def test_display(self):
        """Simply plot the PSF images for the three DUs.
        """
        for du_id in DU_IDS:
            file_path = irf_file_path(self._IRF_NAME, du_id, 'psf')
            psf = xPointSpreadFunction2d(file_path)
            plt.figure('PSF %s DU %d' % (self._IRF_NAME, du_id))
            psf.plot(stretch='log', zlabel='PSF', vmin=1.e-7, vmax=1.e-2)

    @unittest.skip('Waiting for the PSF files.')
    def test_rvs(self, num_events=1000000, du_id=3):
        """Throw random ra and dec offsets for the PSF of one of the DUs.
        """
        file_path = irf_file_path(self._IRF_NAME, du_id, 'psf')
        psf = xPointSpreadFunction2d(file_path)
        ra, dec = psf.delta(num_events)
        plt.figure('PSF %s DU %d rvs' % (self._IRF_NAME, du_id))
        binning = numpy.linspace(-0.1, 0.1, 500)
        hist = xHistogram2d(binning, binning)
        hist.fill(ra, dec)
        hist.plot(logz=True)
        plt.gca().set_aspect('equal')

    @unittest.skip('Waiting for the PSF files.')
    def test_energy_dependence(self, xmax=75.):
        """
        """
        for du_id in DU_IDS:
            plt.figure('PSF energy dependence for DU %d' % du_id)
            for energy in ('229', '451', '640'):
                file_name = 'ixpe_m%d_20210103_psfimage_01_PSF_%s_0_000.fits' % (du_id, energy)
                file_path = os.path.join(IXPEOBSSIM_IRFGEN_DATA, 'mma', file_name)
                psf = xPointSpreadFunction2d(file_path)
                eef, hew = psf.build_eef()
                label = 'PSF 2d @ %.2f keV' % (0.01 * float(energy))
                eef.plot(label=label)
            psf1d = load_psf('ixpe:obssim:v11', du_id)
            plt.plot(psf1d.eef.x, psf1d.eef.y, ls='dashed', label='PSF 1d')
            setup_gca(legend=True, xlabel='r [arcsec]', ylabel='EEF', xmax=xmax, grids=True)

    @unittest.skip('Waiting for the PSF files.')
    def test_radial_profile(self, num_events=1000000):
        """
        """
        for du_id in DU_IDS:
            psf1d = load_psf('ixpe:obssim:v11', du_id)
            file_path = irf_file_path(self._IRF_NAME, du_id, 'psf')
            psf2d = xPointSpreadFunction2d(file_path)
            ra1, dec1 = psf1d.delta(num_events)
            r1 = degrees_to_arcmin(numpy.sqrt(ra1**2. + dec1**2.))
            ra2, dec2 = psf2d.delta(num_events)
            r2 = degrees_to_arcmin(numpy.sqrt(ra2**2. + dec2**2.))
            dx = numpy.random.normal(0., 0.1, size=num_events)
            dy = numpy.random.normal(0., 0.1, size=num_events)
            dr = numpy.sqrt(dx**2. + dy**2.)
            dr = degrees_to_arcmin(numpy.degrees(dr / FOCAL_LENGTH))
            r3 = numpy.sqrt(r2**2. + dr**2.)
            binning = numpy.linspace(0., 8., 250)
            h1 = xHistogram1d(binning).fill(r1)
            h2 = xHistogram1d(binning).fill(r2)
            h3 = xHistogram1d(binning).fill(r3)
            plt.figure('Radial profile DU %d' % du_id)
            h1.plot(label='1d')
            h2.plot(label='2d')
            h3.plot(label='2d + GPD')
            setup_gca(xlabel='r [arcmin]', ylabel='dN/dr', legend=True, logy=True)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

#!/usr/bin/env python
#
# Copyright (C) 2019, the ixpeobssim team.
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
import os
import sys

import numpy
from astropy import wcs

from ixpeobssim.srcmodel.polarization import xStokesSkyMap, xStokesSkyCube
from ixpeobssim import IXPEOBSSIM_CONFIG_FITS
from ixpeobssim.srcmodel.img import xFITSImage
from ixpeobssim.utils.units_ import arcmin_to_degrees
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()


@unittest.skip('Initialization from xy deprecated')
class testStokesSkyCube(unittest.TestCase):

    """Unit test for the Stokes sky cubes.
    """

    def setUp(cls):
        """Create a Stokes sky cube.
        """
        cls.cube = xStokesSkyCube()
        inputs = [
            ('polx_0.4_pf_0.30_radial.fits', 'poly_0.4_pf_0.30_radial.fits', 1., None),
            ('polx_0.4_pf_0.30_radial.fits', 'poly_0.4_pf_0.30_radial.fits', 2., 2.83),
            ('polx_0.4_pf_0.60_radial.fits', 'poly_0.4_pf_0.60_radial.fits', 2.83, 4.),
            ('polx_0.4_pf_0.85_radial.fits', 'poly_0.4_pf_0.85_radial.fits', 4.0, 5.66),
            ('polx_0.4_pf_0.90_radial.fits', 'poly_0.4_pf_0.90_radial.fits', 5.66, 8.),
            ('polx_0.4_pf_0.90_radial.fits', 'poly_0.4_pf_0.90_radial.fits', 12., None)
        ]
        for x_file_name, y_file_name, emin, emax in inputs:
            x_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, x_file_name)
            y_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, y_file_name)
            cls.cube.add_layer_xy(x_file_path, y_file_path, emin, emax)

    def test_cube(self):
        """
        """
        self.cube.plot()
        ra0, dec0 = self.cube.center
        print(ra0, dec0)
        energy = numpy.linspace(1., 12., 100)

        plt.figure('Polarization degree')
        for delta in [1., 2., 3., 4.]:
            ra = numpy.full(energy.shape, ra0)
            dec = numpy.full(energy.shape, dec0 + arcmin_to_degrees(delta))
            pd = self.cube.polarization_degree(ra, dec, energy)
            plt.plot(energy, pd, label='@ (%.4f, %.4f)' % (ra[0], dec[0]))
        plt.xlabel('Energy [keV]')
        plt.ylabel('Polarization degree')
        plt.grid(which='both')
        plt.legend()

    def test_clipping(self, num_events=100000):
        """
        """
        ra0, dec0 = self.cube.center
        ra = numpy.full(num_events, ra0)
        dec = numpy.full(num_events, dec0)
        delta = numpy.random.uniform(-8., 8., num_events)
        dec += arcmin_to_degrees(delta)
        energy = numpy.random.uniform(0., 20., num_events)
        q, u = self.cube(ra, dec, energy)

    def test_tycho(self):
        """
        """
        file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'tycho_4p1_6p1_keV.fits')
        plt.figure('Test Tycho arrows')
        img = xFITSImage(file_path)
        fig = img.plot()
        x_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'polx_0.4_pf_0.90_radial.fits')
        y_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'poly_0.4_pf_0.90_radial.fits')
        pol_map = xStokesSkyMap.load_from_xy(x_file_path, y_file_path)
        pol_map.plot_arrows(nside=40, threshold=0.001)

    def test_pad(self):
        """
        """
        cube = xStokesSkyCube()
        inputs = [
            ('MSH_Rad_PF.fits', 'MSH_Rad_PA.fits', 1., None),
            ('MSH_Rad_PF.fits', 'MSH_Rad_PA.fits', 12., None)
        ]
        for pd_file_name, pa_file_name, emin, emax in inputs:
            pd_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, pd_file_name)
            pa_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, pa_file_name)
            cube.add_layer_pda(pd_file_path, pa_file_path, emin, emax)
        layer = cube.layer(0)
        file_name = 'msh1552_chandra_E2100-10000_FLUXED.fits'
        file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, file_name)
        img = xFITSImage(file_path)
        plt.figure('MSH 1552')
        fig = img.plot()
        layer.plot_arrows(radius=5.)
        plt.figure('MSH 1552 polarization degree')
        layer.plot_polarization_degree()
        layer.plot_arrows(scale=7.5, radius=5.)
        plt.figure('MSH 1552 polarization angle')
        layer.plot_polarization_angle()



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

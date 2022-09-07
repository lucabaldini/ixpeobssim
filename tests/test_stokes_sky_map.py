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

from ixpeobssim.srcmodel.polarization import xStokesSkyMap
from ixpeobssim import IXPEOBSSIM_CONFIG
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.srcmodel.img import xFITSImage
from ixpeobssim.core.stokes import xModelStokesParameters

if sys.flags.interactive:
    plt.ion()



class testStokesSkyMap(unittest.TestCase):

    """Unit test for the Stokes sky maps.
    """

    def setUp(cls):
        """Setup a dummy wcs.
        """
        cls.ra = 350.866590103
        cls.dec = 58.8117384439
        cls.nside = 50
        cls.wcs = wcs.WCS(naxis=2)
        cls.wcs.wcs.crpix = [25., 25.]
        cls.wcs.wcs.cdelt = numpy.array([-0.00279893338726, 0.00279893338726])
        cls.wcs.wcs.crval = [cls.ra, cls.dec]
        cls.wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        cls.wcs.pixel_shape = (cls.nside, cls.nside)

    def test_constant(self, q=0.5, u=0.5):
        """Test a constant Stokes sky map built programmatically from the
        constructor.
        """
        qdata = numpy.full((self.nside, self.nside), q)
        udata = numpy.full((self.nside, self.nside), u)
        map_ = xStokesSkyMap(qdata, udata, self.wcs, None, None)
        _q , _u = map_(self.ra, self.dec)
        self.assertAlmostEqual(_q, q)
        self.assertAlmostEqual(_u, u)
        pd = map_.polarization_degree(self.ra, self.dec)
        pa = map_.polarization_angle(self.ra, self.dec)
        self.assertAlmostEqual(pd, xModelStokesParameters.polarization_degree(q, u))
        self.assertAlmostEqual(pa, xModelStokesParameters.polarization_angle(q, u))

    @unittest.skip('Initialization from xy deprecated')
    def test_xy_tycho(self):
        """This is the unit tests that we used for debugging the Stokes
        sky-map machinery.
        """
        # Paths to the relevant files.
        folder = os.path.join(IXPEOBSSIM_CONFIG, 'fits')
        file_path_x = os.path.join(folder, 'polx_0.4_pf_0.90_radial.fits')
        file_path_y = os.path.join(folder, 'poly_0.4_pf_0.90_radial.fits')
        # Read the input files with the polarization components as images.
        img_x = xFITSImage(file_path_x)
        img_y = xFITSImage(file_path_y)
        # Create the glorious Stokes sky-map.
        map_ = xStokesSkyMap.load_from_xy(file_path_x, file_path_y)

        for i, j in [(15, 29), (28, 15), (24, 23)]:
            ra, dec = img_x.wcs.wcs_pix2world(i, j, True)
            # Read the components from the imput images.
            px = img_x(i, j)
            py = img_y(i, j)
            # Swap, because the input maps are rotated.
            py, px = px, -py
            # Calculate all the input data
            q, u = xModelStokesParameters.xy_to_qu(px, py)
            pd_mod = xModelStokesParameters.polarization_degree(q, u)
            pa_mod = xModelStokesParameters.polarization_angle(q, u)
            print('(i, j) = (%d, %d), (Ra, Dec) = (%.5f, %.5f)' % (i, j, ra, dec))
            print('---Input data:')
            print('px = %.5f, py = %.5f' % (px, py))
            print('q = %.5f, u = %.5f' % (q, u))
            print('pd = %.5f, pa = %.2f deg' % (pd_mod, numpy.degrees(pa_mod)))
            # And now interpolate.
            q, u = map_(ra, dec)
            pd_int = xModelStokesParameters.polarization_degree(q, u)
            pa_int = xModelStokesParameters.polarization_angle(q, u)
            print('---Interpolated data:')
            print('q = %.5f, u = %.5f' % (q, u))
            print('pd = %.5f, pa = %.2f deg' % (pd_int, numpy.degrees(pa_int)))
            delta = abs(pd_mod - pd_int) / pd_mod
            self.assertTrue(delta <= 0.05)

        map_.plot_input_data()
        plt.figure('Polarization degree')
        map_.plot_polarization_degree()
        map_.plot_arrows(nside=50)
        plt.figure('Polarization angle')
        map_.plot_polarization_angle()
        plt.figure('Stokes Q')
        map_.plot_q()
        plt.figure('Stokes U')
        map_.plot_u()



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

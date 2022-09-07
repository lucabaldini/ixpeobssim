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

import unittest
import sys
import os

import numpy

from ixpeobssim import IXPEOBSSIM_CONFIG_FITS, IXPEOBSSIM_TEST_DATA
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.core.fitsio import xFITSImageBase
from ixpeobssim.core.fitsio import xBinTableHDUBase, xPrimaryHDU


if sys.flags.interactive:
    plt.ion()


class TestFitsio(unittest.TestCase):

    """Unit test for the core.fitsio module.
    """

    @classmethod
    def setUpClass(cls):
        """Load an underlying image.
        """
        file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'crab_0p3_10p0_keV.fits')
        cls.image = xFITSImageBase(file_path)

    def test_misc(self):
        """Test the native binary table.
        """
        class CustomBinTable(xBinTableHDUBase):

            NAME = 'CUSTOM'
            HEADER_KEYWORDS = [
                ('PKEY', None)
            ]
            HEADER_COMMENTS = [
                'A comment'
            ]
            DATA_SPECS = [
                ('ENERGY' , 'E', 'keV'  , 'The energy value'),
                ('CHANNEL', 'I'),
                ('AEFF'   , 'E', 'cm**2')
            ]

        hdu = xPrimaryHDU()
        print(hdu)
        n = 20
        data = [numpy.linspace(1., 10., n), numpy.arange(n), numpy.ones(n)]
        table1 = CustomBinTable()
        print(table1)
        keywords = [
            ('AKEY1', 0, 'A keyword'),
            ('AKEY2', 1)
        ]
        comments = ['Howdy, partner?']
        table2 = CustomBinTable(data, keywords, comments)
        print(table2)

    def test_plot(self):
        """Test plotting.
        """
        plt.figure('Crab image test')
        self.image.plot()

    def test_recenter(self, radius=45.):
        """
        """
        plt.figure('Crab image recentered')
        self.image.plot()
        self.image.recenter(83.633083, 22.014500, radius)
        self.image.add_circle(83.633083, 22.014500, radius)
        print(self.image.sky_bounding_box())

    def test_vset(self):
        """Test setting vmin and vmax.
        """
        plt.figure('Crab image adjusted')
        self.image.plot(vmin=-10., vmax=50.)

    def test_stretch(self):
        """Test stretches.
        """
        for stretch in ['log', 'sqrt']:
            plt.figure('Crab image %s' % stretch)
            self.image.plot(stretch=stretch)

    def test_label(self):
        """Test displaying an image with a label.
        """
        plt.figure('Crab image label')
        self.image.plot()
        self.image.add_label('Crab')

    def test_coordinates(self):
        """
        """
        file_path = os.path.join(IXPEOBSSIM_TEST_DATA, 'point_count_map.fits')
        image = xFITSImageBase(file_path)
        ra0, dec0 = image.center()
        dra, ddec = image.delta()
        # Check the image shape.
        self.assertEqual(image.shape(), (200, 200))
        # Make sure the center pixel has signal.
        self.assertTrue(image(100, 100) > 0)
        # The test image is centered at (45., 45.)
        self.assertAlmostEqual(ra0, 45.)
        self.assertAlmostEqual(dec0, 45.)
        self.assertAlmostEqual(image.inner_radius(), 0.027777777777778)

    def test_arrows(self):
        """
        """
        file_path = os.path.join(IXPEOBSSIM_TEST_DATA, 'point_count_map.fits')
        image = xFITSImageBase(file_path)
        plt.figure('Point source with arrows')
        image.plot()
        model = lambda x, y : (numpy.full(x.shape, 1.), numpy.full(y.shape, 1.))
        image.plot_arrows(20, model)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

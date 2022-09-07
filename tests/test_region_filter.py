
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

import unittest
import sys

import ixpeobssim.utils.astro

import ixpeobssim.utils.astro  as astro

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy
import regions

from ixpeobssim.core.fitsio import xFITSImageBase
from ixpeobssim.config.toy_gauss_disk import ROI_MODEL as input_model
from ixpeobssim.core.hist import xHistogram2d
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.evt.fmt import _SKYCOORD_NUM_SIDE_PIXELS
from ixpeobssim.utils.astro import ds9_region_filter_xy, ds9_region_filter_sky, wcs_digitize
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()


class TestRegionFilter(unittest.TestCase):

    """Unit test for the new astropy region filtering.
    """

    @classmethod
    def setUpClass(cls):
        """Run a small simulation with a toy disk and grab the event file.
        """
        pipeline.reset('toy_gauss_disk', overwrite=True)
        file_list = pipeline.xpobssim(duration=100000.)
        file_path = file_list[0]
        cls.event_file = xEventFile(file_path)

    def test_xy(self):
        """Test the filtering in pixel coorinates.
        """
        x, y = self.event_file.xy_data()
        pix = 0.5 * (_SKYCOORD_NUM_SIDE_PIXELS - 1)
        center = regions.PixCoord(pix, pix)
        circle = regions.CirclePixelRegion(center, 50.)
        rectangle = regions.RectanglePixelRegion(center=center, width=50., height=200.)
        region_list = [circle, rectangle]

        binning = numpy.linspace(0, _SKYCOORD_NUM_SIDE_PIXELS - 1, _SKYCOORD_NUM_SIDE_PIXELS)

        def plot_hist(mask, figure_name):
            """Small nested function to plot xy maps with various filters.
            """
            plt.figure('%s (xy)' % figure_name)
            if mask is None:
                mask = numpy.ones(x.shape, dtype=bool)
            logger.info('xy %s -> %d event(s)', figure_name, mask.sum())
            xHistogram2d(binning, binning, xlabel='x', ylabel='y').fill(x[mask], y[mask]).plot()

        for mask, figure_name in (
            (None, 'Full'),
            (ds9_region_filter_xy(x, y, circle), 'Circle'),
            (ds9_region_filter_xy(x, y, *region_list), 'Logical OR'),
            (ds9_region_filter_xy(x, y, *region_list, invert=True), 'Logical OR inverted'),
            (ds9_region_filter_xy(x, y, *region_list, compound_mode='and'), 'Logical AND'),
            (ds9_region_filter_xy(x, y, *region_list, compound_mode='xor'), 'Logical XOR'),):
            plot_hist(mask, figure_name)

    def test_sky(self):
        """Test the filtering in sky coordinates.
        """
        wcs_ = self.event_file._wcs
        ra, dec = self.event_file.sky_position_data()
        center = SkyCoord(input_model.ra * units.deg, input_model.dec * units.deg)
        circle = regions.CircleSkyRegion(center, 50. * units.arcsec)
        rectangle = regions.RectangleSkyRegion(center=center, width=50. * units.arcsec, height=200. * units.arcsec)
        region_list = [circle, rectangle]

        def plot_skymap(mask, figure_name):
            """
            """
            plt.figure('%s (sky)' % figure_name)
            if mask is None:
                mask = numpy.ones(ra.shape, dtype=bool)
            logger.info('sky %s -> %d event(s)', figure_name, mask.sum())
            data = wcs_digitize(wcs_, ra[mask], dec[mask])
            xFITSImageBase.make_plot(data, wcs_)

        for mask, figure_name in (
            (None, 'Full'),
            (ds9_region_filter_sky(ra, dec, wcs_, circle), 'Circle'),
            (ds9_region_filter_sky(ra, dec, wcs_, *region_list), 'Logical OR'),
            (ds9_region_filter_sky(ra, dec, wcs_, *region_list, invert=True), 'Logical OR inverted'),
            (ds9_region_filter_sky(ra, dec, wcs_, *region_list, compound_mode='and'), 'Logical AND'),
            (ds9_region_filter_sky(ra, dec, wcs_, *region_list, compound_mode='xor'), 'Logical XOR'),):
            plot_skymap(mask, figure_name)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

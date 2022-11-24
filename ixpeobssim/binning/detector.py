# Copyright (C) 2015--2022, the ixpeobssim team.
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

"""Binning data products in detector coordinates.
"""

from __future__ import print_function, division

from astropy.io import fits
import numpy

from ixpeobssim.binning.base import xEventBinningBase, xBinnedFileBase
from ixpeobssim.core.hist import xGpdMap2d
from ixpeobssim.instrument.gpd import gpd_map_binning, GPD_PHYSICAL_HALF_SIDE_X,\
    GPD_PHYSICAL_HALF_SIDE_Y
from ixpeobssim.utils.logging_ import logger


# pylint: disable=invalid-name, attribute-defined-outside-init, too-few-public-methods
# pylint: disable=no-member, arguments-differ


class xEventBinningARMAP(xEventBinningBase):

    """Class for ARMAP binning.
    """

    INTENT = 'rate per unit area map in detector coordinates'
    SUPPORTED_KWARGS = ['npix']

    def process_data(self):
        """Convenience function factoring out the code in common with the
        corresponding EFLUX class---see the overloaded method in there.

        Here we are binning the event position in detector coordinates
        and dividing by the livetime and the bin area. The function returns
        a n x n array of area rate values to be written in the output file.
        """
        detx, dety = self.event_file.det_position_data()
        xbinning, ybinning = gpd_map_binning(GPD_PHYSICAL_HALF_SIDE_X,
            GPD_PHYSICAL_HALF_SIDE_Y, self.get('npix'))
        bin_area = (xbinning[1] - xbinning[0]) * (ybinning[1] - ybinning[0])
        rate, _, _ = numpy.histogram2d(detx, dety, bins=(xbinning, ybinning))
        rate /= self.event_file.livetime() * bin_area
        return rate

    def bin_(self):
        """Overloaded method.
        """
        rate = self.process_data()
        primary_hdu = self.build_primary_hdu(rate)
        primary_hdu.add_keyword('TOTCNTS', self.event_file.num_events(),
                                'total counts in the original photon list')
        hdu_list = fits.HDUList([primary_hdu])
        self.write_output_file(hdu_list)



class xBinnedAreaRateMap(xBinnedFileBase):

    """Display interface to binned ARMAP files.
    """

    Z_TITLE = 'Dead-time corrected rate [Hz mm$^{-2}$]'

    def _read_data(self):
        """Overloaded method.
        """
        self.data = self.hdu_list['PRIMARY'].data.T
        self.counts = self.hdu_list['PRIMARY'].header['TOTCNTS']

    def __iadd__(self, other):
        """Overloaded method for ARMAP binned data addition.

        Note that, given the peculiarity of this binned data product
        (we are typically interested in the area rate per detector)
        we are doing a weighted average of the input data based on the
        number of counts for each DU in the original event list.
        """
        self._check_iadd(other)
        self.data = (self.data * self.counts + other.data * other.counts) /\
                    (self.counts + other.counts)
        self.counts += other.counts
        return self

    def plot(self):
        """Plot the data.
        """
        threshold = 0.1
        sel = self.data[self.data > threshold * self.data.max()]
        logger.info('Average = %.3f (threshold = %.3f), maximum = %.3f',
                    sel.mean(), threshold, sel.max())
        npix = self.data.shape[0]
        hist = xGpdMap2d(npix, zlabel=self.Z_TITLE)
        hist.set_content(self.data, self.data / self.data.sum() * self.counts)
        hist.plot()



class xEventBinningEFLUX(xEventBinningARMAP):

    """Class for EFLUX binning.
    """

    INTENT = 'energy-flux map in detector coordinates'

    def process_data(self):
        """Overloaded method.

        Here we are just multiplying by the average measured event energy.
        """
        rate = xEventBinningARMAP.process_data(self)
        rate *= self.event_file.energy_data(mc=False).mean()
        return rate



class xBinnedAreaEnergyFluxMap(xBinnedAreaRateMap):

    """Display interface to binned EFLUX files.
    """

    Z_TITLE = 'Dead-time corrected energy flux [keV mm$^{-2}$ s$^{-1}$]'

#!/urs/bin/env python
#
# Copyright (C) 2015--2020, the ixpeobssim team.
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

"""FITS image support.
"""

from __future__ import print_function, division


import numpy

from ixpeobssim.core.fitsio import xFITSImageBase


# pylint: disable=invalid-name, too-many-arguments


class xFITSImage(xFITSImageBase):

    """Class describing a FITS image equipped to extract random coordinates.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xFITSImageBase.__init__(self, file_path)
        self.cdf = self._build_cdf()

    def _build_cdf(self):
        """Build the cumulative distribution function.

        (This is used to extract random positions from the image when
        simulating extended sources.)
        """
        # Note the cast to float is a terrible workaround for issue
        # https://bitbucket.org/ixpesw/ixpeobssim/issues/608
        # triggered by numpy 1.22.0
        cdf = numpy.cumsum(self.data.ravel().astype(float))
        cdf /= cdf[-1]
        return cdf

    def rvs_coordinates(self, size=1, randomize=True):
        """Generate random coordinates based on the image map.

        Arguments
        ---------
        size : int
            The number of sky coordinates to be generated.

        randomize : bool
            If true, the positions are randomized uniformely within each pixel.
        """
        # Throw a vector of random numbers.
        u = numpy.random.rand(size)
        # Convert u into pixel indices through the underlying cdf.
        pixel = numpy.searchsorted(self.cdf, u)
        # Convert from pixel serial id to bi-dimensional coordinates.
        row, col = numpy.unravel_index(pixel, self.data.shape)
        # Stack the array vertically and transpose the result---we need a an
        # N x 2 array, here---and remember that the transpose operator is no-op
        # for one-dimensional numpy arrays.
        pixel_coords = numpy.vstack((col, row)).transpose()
        # Switch to world coordinates through the WCS.
        # Note we are calling wcs_world2pix with origin=0, as the input is read
        # in numpy space---this has been tested with a single point source at
        # the center of the field of view. See
        # https://docs.astropy.org/en/stable/api/astropy.wcs.WCS.html#astropy.wcs.WCS.wcs_world2pix
        world_coords = self.wcs.wcs_pix2world(pixel_coords, 0)
        # Unpack the output into the ra and dec arrays.
        ra, dec = world_coords[:, 0], world_coords[:, 1]
        # If needed, add some randomization.
        if randomize:
            delta_ra = 0.5 * self.primary_hdu.header['CDELT1']
            delta_dec = 0.5 * self.primary_hdu.header['CDELT2']
            ra += numpy.random.uniform(-delta_ra, delta_ra, size)
            dec += numpy.random.uniform(-delta_dec, delta_dec, size)
        return ra, dec

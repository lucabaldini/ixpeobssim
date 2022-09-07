# Copyright (C) 2015--2019, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public License as published by
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

"""PSF parametrization.
"""

from __future__ import print_function, division

from ixpeobssim.core.spline import xInterpolatedBivariateSplineLinear
from ixpeobssim.irf.base import xResponseBase
from ixpeobssim.utils.matplotlib_ import plt

# pylint: disable=invalid-name, no-member


class xVignetting(xResponseBase, xInterpolatedBivariateSplineLinear):

    """Class describing vignetting.

    The vignetting is essentially a bivariate linear spline, with built-in
    facilities for evaluation and plotting.

    Arguments
    ---------
    file_path : str
        The path to the vign.fits file containing the vignettting table.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xResponseBase.__init__(self, file_path, 'fits')
        x = 0.5 * (self.field('ENERG_LO') + self.field('ENERG_HI'))[0]
        y = self.field('THETA')[0]
        z = self.field('VIGNETTING')[0].T
        fmt = dict(xlabel='Energy [keV]', ylabel='Off-axis angle [arcmin]',
            zlabel='Vignetting coefficient')
        xInterpolatedBivariateSplineLinear.__init__(self, x, y, z, **fmt)

    def plot(self):
        """Plot the vignetting.
        """
        #pylint: disable=arguments-differ
        plt.figure(self.base_name)
        xInterpolatedBivariateSplineLinear.plot(self)
        xInterpolatedBivariateSplineLinear.plot_contours(self)

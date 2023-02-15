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

"""PSF parametrization.
"""

from __future__ import print_function, division

from astropy.io import fits
import numpy
from scipy.ndimage.filters import gaussian_filter

from ixpeobssim.irf.base import xResponseBase
from ixpeobssim.core.fitsio import xFITSImageBase
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.core.rand import xUnivariateGenerator
from ixpeobssim.srcmodel.img import xFITSImage
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, labeled_marker
from ixpeobssim.utils.units_ import arcsec_to_degrees, degrees_to_arcsec

# pylint: disable=invalid-name, too-many-arguments, no-member


def gauss_king(r, W, sigma, N, r_c, eta):
    """Functional representation of the Gaussian plus King PSF profile
    described in `Fabiani et al., 2014 <http://arxiv.org/abs/1403.7200>`_,
    equation (2):

    .. math::
        \\text{PSF}(r) = W \\exp^{-(\\frac{r^2}{2\\sigma^2})} +
        N\\left( 1 + \\left( \\frac{r}{r_c} \\right)^2 \\right)^{-\\eta}

    Arguments
    ---------
    r : float or array
        The radial distance from the true source position is arcsec.

    W : float
        Normalization of the Gaussian component.

    sigma : float
        Width of the Gaussian component.

    N : float
        Normalization of the King component.

    r_c : float
        Characteristic radius of the King component.

    eta : float
        Exponent of the King component.
    """
    return W * numpy.exp(-(r**2. / (2. * sigma**2.))) + N * (1. + (r / r_c)**2.)**(-eta)


def gauss_king_eef_at_infinity(W, sigma, N, r_c, eta):
    """Return the value of the Encircled Energy Fraction (EEF) at infinity,
    given the parameters of the functional representation, see equation (4)
    of `Fabiani et al., 2014 <http://arxiv.org/abs/1403.7200>`_.

    .. math::
        \\text{EEF}(\\infty) = 2\\pi W\\sigma^2 +
        \\pi\\frac{r_c^2 N}{\\eta - 1}

    Arguments
    ---------
    r : float or array
        The radial distance from the true source position is arcsec.

    W : float
        Normalization of the Gaussian component.

    sigma : float
        Width of the Gaussian component.

    N : float
        Normalization of the King component.

    r_c : float
        Characteristic radius of the King component.

    eta : float
        Exponent of the King component.
    """
    return 2. * numpy.pi * W * sigma**2. + numpy.pi * r_c**2. * N / (eta - 1.)



class xPointSpreadFunctionBase(xResponseBase):

    """Base virtual class for the PSF data structures.

    .. versionadded:: 28.1.0

       This was added in version 27.1.0 as a base class for the old-style PSF
       (azimuthally symmetric and including the position resolution of the GPD)
       and the new-style PSF necessary to study the Stokes cross-talk
       (non azimuthally symmetric and limited to the X-ray optics).
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xResponseBase.__init__(self, file_path, 'fits')

    def delta(self, size=1):
        """Return an array of random offset (ra, dec) due to the PSF.
        """
        raise NotImplementedError

    def smear(self, ra, dec):
        """Smear a pair of arrays of coordinates.
        """
        assert ra.size == dec.size
        delta_ra, delta_dec = self.delta(ra.size)
        return ra + delta_ra / numpy.cos(numpy.radians(dec)), dec + delta_dec

    def smear_single(self, ra, dec, size=1):
        """Smear a single pair of coordinates for an arbitrary number of times.
        """
        return self.smear(numpy.full(size, ra), numpy.full(size, dec))



class xPointSpreadFunction(xInterpolatedUnivariateSpline, xPointSpreadFunctionBase):

    """Class describing a (simplified, energy independent) PSF.

    The effective area is essentially a linear spline, with built-in facilities
    for evaluation and plotting.

    Arguments
    ---------
    psf_file_path : str
        The path to the .psf FITS file containing the PSF parameters.

    Note
    ----
    The parametrization is taken from `Fabiani et al., 2014
    <http://arxiv.org/abs/1403.7200>`_, table 2.
    """

    MAX_RADIUS = 480.
    PARAM_NAMES = ['W', 'sigma', 'N', 'r_c', 'eta']

    def __init__(self, file_path):
        """Constructor.
        """
        xPointSpreadFunctionBase.__init__(self, file_path)
        data = self.hdu_list['PSF'].data
        W = data['W']
        sigma = data['SIGMA']
        N = data['N']
        r_c = data['R_C']
        eta = data['ETA']
        self.__params = (W, sigma, N, r_c, eta)
        # Tabulate the actual PSF values.
        r = numpy.linspace(0, self.MAX_RADIUS, int(self.MAX_RADIUS) + 1)
        pdf = gauss_king(r, *self.__params)
        fmt = dict(xlabel='r [arcsec]', ylabel='PSF [sr$^{-1}$]')
        xInterpolatedUnivariateSpline.__init__(self, r, pdf, k=2, **fmt)
        # Include the solid angle for the actual underlying random generator.
        pdf *= 2 * numpy.pi * r
        fmt = dict(xlabel='r [arcsec]', ylabel='$2 \\pi r \\times$ PSF')
        self.generator = xUnivariateGenerator(r, pdf, k=1, **fmt)
        # Finally, calculate the EEF and HEW.
        self.eef, self.hew = self.build_eef()

    def build_eef(self):
        """Build the Encircled Energy Fraction (EEF) as a function of r.

        And, while we're at it, we also calculate and cache the HEW.
        """
        _r = self.x
        _y = numpy.array([self.generator.integral(_r[0], _rp) for _rp in _r])
        _y /= gauss_king_eef_at_infinity(*self.__params)
        hew = 2 * xInterpolatedUnivariateSpline(_y, _r, k=2)(0.5)
        fmt = dict(xlabel='r [arcsec]', ylabel='PSF EEF')
        return xInterpolatedUnivariateSpline(_r, _y, k=1, **fmt), hew

    def plot(self):
        """Overloaded plot method.
        """
        # pylint: disable=arguments-differ, unexpected-keyword-arg
        plt.figure('%s PSF' % self.base_name)
        xInterpolatedUnivariateSpline.plot(self, logy=True, label=self.base_name)
        plt.axis([0, self.MAX_RADIUS, None, None])
        plt.grid(which='both')
        plt.legend()
        plt.figure('%s PSF EEF' % self.base_name)
        self.eef.plot(label=self.base_name)
        x = self.hew / 2.
        y = self.eef(x)
        label = 'HEW = %.2f arcsec' % self.hew
        labeled_marker(x, y, label, dx=2., va='top')
        plt.axis([0, self.MAX_RADIUS, None, None])
        plt.grid(which='both')
        plt.legend()

    def delta(self, size=1):
        """Return an array of random offset (in ra, dec or L, B) due to the PSF.

        Note the output is converted in degrees.
        """
        r = arcsec_to_degrees(self.generator.rvs(size))
        phi = numpy.random.uniform(0., 2. * numpy.pi, size)
        return r * numpy.cos(phi), r * numpy.sin(phi)

    def draw_psf_circle(self, image, x=0.8, y=0.8, label='HPD', hpd_value=False,
        color='white', lw=1.5):
        """Add the PSF circle to the image with labels.

        Arguments
        ---------
        image : xFITSImageBase instance
            The parent xFITSImageBase object

        x, y : float
            The position of the psf circle in relative canvas coordinates
        """
        if hpd_value:
            label = '%s (%.1f")' % (label, self.hew)
        image.add_label(label, x, y + 0.05, color=color, ha='center')
        image.add_circle(x, y, 0.5 * self.hew, mode='canvas', color=color, lw=lw)



class xPointSpreadFunction2d(xPointSpreadFunctionBase):

    """Two-dimensional version (i.e., non azimuthally symmetric) version of the PSF.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xPointSpreadFunctionBase.__init__(self, file_path)
        img_hdu = self.hdu_list['IMG_2D']
        self.psf_image = xFITSImage(img_hdu)

    def delta(self, size=1):
        """Return an array of random offset (in ra, dec) due to the PSF.
        """
        delta_ra, delta_dec = self.psf_image.rvs_coordinates(size)
        # Note that, since the WCS is centered in (0, 0) negative ra values
        # get carried away near 360 degrees, so we have to fold them back
        # explicitly.
        delta_ra[delta_ra > 180.] -= 360.
        return delta_ra, delta_dec



class xPointSpreadFunction4d(xPointSpreadFunction2d):

    """Glorious, 4-dimensional PSF containing an approximate scaling of the
    radius with energy and off-axis angle.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xPointSpreadFunction2d.__init__(self, file_path)

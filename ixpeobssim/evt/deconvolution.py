#!/usr/bin/env python
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

"""Module encapsulating image deconvolution facilities.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.irf import DEFAULT_IRF_NAME, load_psf
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.units_ import degrees_to_arcmin


# pylint: disable=invalid-name


def circular_kernel(nside, normalize=False):
    """Build a simple circular kernel of a given side.

    Arguments
    ---------
    nside : int
        The side of the kernel array, assumed to be squared.

    normalize : bool
        If True, normalize to the kernel content (i.e., use for averaging rather
        than summing).
    """
    shape = (nside, nside)
    kernel = numpy.zeros(shape, dtype=float)
    x, y = numpy.indices(shape) - nside / 2. + 0.5
    kernel[x**2.  + y**2. <= nside**2. / 4.] = 1.
    if normalize:
        kernel /= kernel.sum()
    return kernel


def calculate_psf_kernel(pixel_size, irf_name=DEFAULT_IRF_NAME, du_id=1,
    max_radius=0.03, sample_size=5000000):
    """Calculate the digital PSF kernel for image convolution.

    This is essentially creating a square grid whose pixel size is set to the
    value of the corresponding function argument and whose half size is (loosely)
    determined by the max_radius arguments. A number of points is thrown in the
    grid according to the PSF shape and the corresponding probability for each
    pixel is estimated by numerical integration.

    Arguments
    ---------
    pixel_size : float
        The pixels size of the image that the kernel is operating on
        (in decimal degrees).

    irf_name : str
        The IRF name.

    du_id : int
        The DU id.

    max_radius : float
        The maximum radius (in decimal degrees) for defining the size of the
        kernel grid. This should be large enough to fully contain the PSF.

    sample_size : int
        The size of the random sample for the Monte Carlo integration of the
        psf profile over the square grid.
    """
    logger.info('Building the PSF kernel for "%s" (DU %d)', irf_name, du_id)
    # Note the number of bins on both axes must be odd, so that the origin
    # of the PSF is at the center of a bin.
    nside = 2 * int(max_radius / pixel_size) + 1
    delta = 0.5 * pixel_size * nside
    binning = numpy.linspace(-delta, delta, nside + 1)
    delta = degrees_to_arcmin(delta)
    args = nside, nside, delta
    logger.info('Sampling over a %d x %d spatial grid (+/-%.2f arcmin)', *args)
    # Load the PSF.
    psf = load_psf(irf_name, du_id)
    # Throw a whole nunch of random numbers with the PSF shape and bin them
    # over the kernel.
    delta = psf.delta(sample_size)
    kernel, _ = numpy.histogramdd(delta, (binning, binning))
    # Debug info to keep track of the kernel containment.
    kernel_sum = kernel.sum()
    containment = kernel.sum() / sample_size
    logger.info('Estimated PSF containment: %.9f', containment)
    # Normalize the kernel and we're done.
    return kernel / kernel_sum


def _kernel_radius(kernel):
    """Return the kernel radius in pixels.

    This is meant to return n, where the kernel shape is (2n + 1) x (2n + 1).
    """
    nx, ny = kernel.shape
    assert nx == ny
    return nx // 2


def _kernel_mask(i, j, nk, image_shape):
    """Return the appropriate slice for convolving the generic pixel (i, j)
    within an image. For (i, j) = (50, 50) and a kernel radius of 2 pixels, e.g.,
    this returns (slice(48, 53, None), slice(48, 53, None)).

    It is assumed that the kernel is square, with dimension (2nk + 1) x (2nk + 1).
    """
    # Calculate the slice for the image---this is (i-nk--i+nk), (j-nk--j+nk)
    # modulo border effects.
    nx, ny = image_shape
    imin = max(i - nk, 0)
    imax = min(i + nk + 1, nx)
    jmin = max(j - nk, 0)
    jmax = min(j + nk + 1, ny)
    img_mask = (slice(imin, imax), slice(jmin, jmax))
    # And now the slice for the kernel---this is (0--2nk+1), (0--2nk+1)
    # modulo border effects.
    nkmax = 2 * nk + 1
    imin = max(0, nk - i)
    imax = min(nkmax, nx - i + nk)
    jmin = max(0, nk - j)
    jmax = min(nkmax, ny - j + nk)
    ker_mask = (slice(imin, imax), slice(jmin, jmax))
    return img_mask, ker_mask


def convolve_image(image, kernel):
    """Convolve an image with the PSF kernel.

    While I am sure there exist a more clever (and faster) way to do this with
    scipy, I retained the explicit Python loop, here, as I don't think this is
    a bottleneck in any way, and having a full implementation is useful for
    debugging.

    Arguments
    ---------
    image : 2d array
        The input image to be convolved with the PSF.

    kernel : 2d array
        The PSF kernel for the convolution.
    """
    output = numpy.full(image.shape, 0.)
    nk = _kernel_radius(kernel)
    for i, j in numpy.ndindex(image.shape):
        imask, kmask = _kernel_mask(i, j, nk, image.shape)
        output[imask] += kernel[kmask] * image[i, j]
    return output


def calculate_inverse_kernel(intensity, true_intensity, kernel):
    """Calculate the gargantuan (nx x ny x nk x nk) matrix that allows to
    deconvolve an (I, U or Q) image, given a proxy of the true intensity.

    The basic idea, here, is that if you have a sensible proxy of the true I
    map (e.g., a high-resolution image from Chandra) you can use it to gauge
    the fraction of the measured intensity ending up in each true pixel, and
    you can actually use that to deconvole a map (with the same binning) of
    either Q or U. Phrased in a slighlty different way, while the PSF kernel
    for the direct true -> measured convolution is identical for all pixels, the
    kernel for the inverse measured -> true deconvolution is in general
    different for all the pixels and depends on the true intensity. This is
    why we need to precompute a store a specific kernel for each pixels. This
    can be in turn use to deconvole any measured quantity with the same spatial
    binning (Q, U, W2).

    Warning
    -------
    I am not positive this is right for the border pixels, where the kernel mask
    ends up out of the image---need to study that.

    Arguments
    ---------
    intensity : 2d array
        The measured intensity map.

    true_intensity : 2d array
        The proxy for the true intensity map.

    kernel : 2d array
        The PSF kernel for the convolution.
    """
    # Before we start, we normalize the total intensity of the true image to that
    # of the reference image, since the two might, in principle, have totally
    # arbitrary normalizations.
    true_intensity *= intensity.sum() / true_intensity.sum()
    nk = _kernel_radius(kernel)
    K = numpy.full((*intensity.shape, *kernel.shape), 0.)
    for i, j in numpy.ndindex(intensity.shape):
        imask, kmask = _kernel_mask(i, j, nk, intensity.shape)
        subimg = intensity[imask]
        # Need to protect against zero-division errors.
        mask = (subimg != 0.)
        K[i, j][kmask][mask] = (kernel[kmask][mask] * true_intensity[i, j] / subimg[mask])
    return K


def deconvolve_image(image, kernel, K):
    """Deconvolve a generic image (e.g., Q or U).

    Arguments
    ---------
    image : 2d array
        The input image to be deconvoled.

    kernel : 2d array
        The PSF kernel.

    K : 4d array
        The inverse pixel-by-pixel kernel for the deconvolution.
    """
    output = numpy.full(image.shape, 0.)
    nk = _kernel_radius(kernel)
    for i, j in numpy.ndindex(image.shape):
        imask, kmask = _kernel_mask(i, j, nk, image.shape)
        output[i, j] = (image[imask] * K[i, j][kmask]).sum()
    return output

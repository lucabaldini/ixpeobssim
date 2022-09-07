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


"""Unit test for the evt.deconvolution module.
"""

import unittest
import sys

import numpy
import matplotlib

from ixpeobssim.evt.deconvolution import _kernel_mask, calculate_psf_kernel,\
    convolve_image, deconvolve_image, calculate_inverse_kernel, circular_kernel
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca

if sys.flags.interactive:
    plt.ion()


class TestDeconvolution(unittest.TestCase):

    """Unit test for the evt.deconvolution module.
    """

    @classmethod
    def setUpClass(cls):
        """
        """
        cls.kernel = calculate_psf_kernel(0.00145833333333)
        plt.figure('Kernel')
        plt.imshow(cls.kernel, norm=matplotlib.colors.LogNorm())
        setup_gca(xlabel='$x$ offset [pixels]', ylabel='$y$ offset [pixels]')
        cbar = plt.colorbar()
        cbar.set_label('Probability')

    def test_kernel_mask(self, nk=1, image_shape=(5, 5)):
        """
        """
        for i, j in numpy.ndindex(image_shape):
            img_mask, ker_mask = _kernel_mask(i, j, nk, image_shape)
            print(i, j)
            print(img_mask)
            print(ker_mask)

    def test_circular_kernel(self, nside=15):
        """
        """
        kernel = circular_kernel(nside)
        plt.matshow(kernel)
        plt.colorbar()

    def test_roundtrip(self):
        """Make sure that the convolution-deconvolution process roundtrips
        on a test image.
        """
        img = numpy.full((200, 200), 0.)
        img[100, 100] = 100.
        img[75, 45:110] = 50
        conv = convolve_image(img, self.kernel)
        K = calculate_inverse_kernel(conv, img, self.kernel)
        deconv = deconvolve_image(conv, self.kernel, K)
        plt.figure('Original image')
        plt.imshow(img)
        plt.colorbar()
        plt.figure('Convolved image')
        plt.imshow(conv)
        plt.colorbar()
        plt.figure('Deconvolved image')
        plt.imshow(deconv)
        plt.colorbar()
        plt.figure('Image difference')
        plt.imshow(deconv - img)
        plt.colorbar()
        self.assertTrue(numpy.allclose(img, deconv))
        max_diff = abs(deconv - img).max()
        logger.info('Maximum image difference is %.3e', max_diff)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

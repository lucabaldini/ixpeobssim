#!/usr/bin/env python
#
# Copyright (C) 2018, the ixpeobssim team.
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

import os
import unittest

from astropy.io import fits
import numpy

from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.utils.environment import NUMPY_VERSION
from ixpeobssim.utils.logging_ import logger


class TestIssue(unittest.TestCase):

    """Unit test for issue 608 https://bitbucket.org/ixpesw/ixpeobssim/issues/608
    """

    @classmethod
    def setUpClass(cls):
        """Open the FITS file for the Cal C morphology definition and retrieve the
        data.
        """
        file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'fits', 'onboard_cal_C_img.fits')
        logger.info('Opening input FITS file %s...', file_path)
        with fits.open(file_path) as hdu_list:
            hdu_list.info()
            cls.data = hdu_list['PRIMARY'].data
            shape = cls.data.shape
            logger.info('Input image data: %s (shape %s)', cls.data, shape)
            logger.info('Data dtype: %s', cls.data.dtype)

    def test_cast(self):
        """Cast the data to float and make sure we're not changing the values.
        """
        self.assertTrue(numpy.allclose(self.data, self.data.astype(float)))
        self.assertTrue(numpy.allclose(self.data.ravel(), self.data.ravel().astype(float)))

    def _calculate_cdf(self):
        """Helper function to verify the numpy exception.
        """
        return numpy.cumsum(self.data.ravel())

    def test_issue_608(self):
        """Small test for issue 608.
        """
        if NUMPY_VERSION == '1.22.0':
            logger.info('numpy version 1.22.0 will trigger issue 608...')
            self.assertRaises(TypeError, self._calculate_cdf)
            cdf = numpy.cumsum(self.data.ravel().astype(float))
        else:
            cdf = numpy.cumsum(self.data.ravel())



if __name__ == '__main__':
    unittest.main()

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


"""Unit test for packaging information.
"""


import unittest

import numpy
import scipy
import astropy
import matplotlib

from ixpeobssim.utils.packaging_ import xPackageVersion, parse_version_string, retrieve_version


class testPackaging(unittest.TestCase):

    """Unit test for the extraction of packaging information.
    """

    def test_main_packages(self):
        """Print the version of our major dependencies.
        """
        for package in (numpy, scipy, astropy, matplotlib):
            print(package.__name__, retrieve_version(package))

    def test_strings(self):
        """Make sure we parse version string without patch and with post suffix.
        """
        print(parse_version_string('1.2'))
        print(parse_version_string('4.0.1.post1'))

    def test_xspec_patch(self):
        """See issue https://bitbucket.org/ixpesw/ixpeobssim/issues/507/
        """
        with self.assertRaises(SystemExit):
            parse_version_string('12.12.0g')

    def test_comparisons(self):
        """
        """
        self.assertTrue(xPackageVersion(1, 4, 0) < '1.6.0')
        self.assertTrue(xPackageVersion(1, 4, 0) <= '1.6.0')
        self.assertTrue(xPackageVersion(1, 4, 0) != '1.6.0')
        self.assertTrue(xPackageVersion(1, 4, 0) != '1.6')
        self.assertTrue(xPackageVersion(1, 4, 0) != '1.6.0.post1')
        self.assertTrue(xPackageVersion(1, 4, 0) == '1.4.0')
        self.assertTrue(xPackageVersion(1, 4, 0) <= '1.4.0')
        self.assertTrue(xPackageVersion(1, 6, 0) >= '1.4.0')
        self.assertTrue(xPackageVersion(1, 6, 0) >= '0.4')
        self.assertTrue(xPackageVersion(1, 6) >= '0.4.0')



if __name__ == '__main__':
    unittest.main()

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

"""Unit tests for Python environment.
"""

from __future__ import print_function, division

import os
import unittest
import re
import importlib

import ixpeobssim
from ixpeobssim.utils.packaging_ import retrieve_version
from ixpeobssim.utils.logging_ import logger


class TestEnvironment(unittest.TestCase):

    """Unit test for the ixpeobssim.utils.system_ module.
    """

    def test_python(self):
        """Test the Python version.
        """
        logger.info('Checking Python...')
        logger.info('%s installed', ixpeobssim.PYTHON_VERSION)
        self.assertTrue(ixpeobssim.PYTHON_VERSION >= '3.6')

    def test_deps(self):
        """Test the dependencies.
        """
        logger.info('Checking dependencies...')
        file_path = os.path.join(ixpeobssim.IXPEOBSSIM_BASE, 'requirements.txt')
        with open(file_path) as input_file:
            for line in input_file:
                line = line.strip('\n')
                try:
                    package, operator, version = re.split('(>=|==|<=|>|<)', line)
                    logger.info('Checking %s...', package)
                    package = importlib.import_module(package)
                    installed = retrieve_version(package)
                    logger.info('%s installed', installed)
                    self.assertTrue(eval('installed %s version' % operator))
                except ValueError:
                    pass


if __name__ == '__main__':
    unittest.main()

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

import glob
import os
import unittest

from ixpeobssim import IXPEOBSSIM_BIN
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.system_ import import_module


class TestMain(unittest.TestCase):

    """Make sure all the apps in the bin folder have the main() entry point
    properly defined. This is necessary to call the python apps without
    explicitely including the file extension when the package is installed
    (i.e., in user mode).

    See https://bitbucket.org/ixpesw/ixpeobssim/issues/480 and
    https://bitbucket.org/ixpesw/ixpeobssim/issues/480
    """

    def test_apps(self):
        """
        """
        for file_path in glob.glob(os.path.join(IXPEOBSSIM_BIN, 'xp*.py')):
            logger.info('Testing main() hook for %s...', file_path)
            module = import_module(file_path)
            msg = '%s has no main hook' % module
            self.assertTrue('main' in dir(module), msg)



if __name__ == '__main__':
    unittest.main()

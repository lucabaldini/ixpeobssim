#!/usr/bin/env python
#
# Copyright (C) 2019, the ixpeobssim team.
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


"""Dummy unit test for the logger
"""


import unittest

from ixpeobssim.utils.logging_ import logger


class testLogger(unittest.TestCase):

    """Unit test for the logger.
    """

    def test(self):
        """
        """
        logger.debug('Debug message...')
        logger.info('Info message...')
        logger.warning('Warning message')
        logger.error('Error message...')
        logger.debug('Debug message...')


if __name__ == '__main__':
    unittest.main()

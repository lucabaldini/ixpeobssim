#!/usr/bin/env python
#
# Copyright (C) 2016, the ixpeobssim team.
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


"""Test the matplotlib settings.
"""


import unittest

import ixpeobssim.utils.os_ as os_

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.utils.matplotlib_ import plt, save_gcf



class testMatplotlib(unittest.TestCase):

    """Unit test for the matplotlib_ module.
    """

    def test_save_gcf(self):
        """
        """
        file_list = []
        plt.figure('Test figure')
        plt.plot([1, 2, 3], [1, 4, 9])
        file_list += save_gcf(IXPEOBSSIM_DATA)
        file_list += save_gcf(IXPEOBSSIM_DATA, 'test')
        file_list += save_gcf(IXPEOBSSIM_DATA, 'test', ['png', 'pdf'])
        for file_path in file_list:
            os_.rm(file_path)



if __name__ == '__main__':
    unittest.main()

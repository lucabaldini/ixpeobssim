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


"""Unit tests for the ixpeobssim.utils.os_ module.
"""


import unittest
import tempfile
import os

from ixpeobssim.utils.os_ import cp, mv, rm, mkdir, rmdir


class testos(unittest.TestCase):

    """Unit test for the ixpeobssim.utils.os_ module.
    """

    def test_cp_mv_rm(self):
        """Test copying, renaming and moving files...
        """
        f = tempfile.NamedTemporaryFile(delete=False)
        f.close()
        src = f.name
        dest1 = '%s_copy' % src
        dest2 = '%s_new' % src
        self.assertEqual(cp(src, dest1), 0)
        self.assertEqual(mv(dest1, dest2), 0)
        self.assertEqual(rm(dest1), 0)
        self.assertEqual(rm(dest2), 0)
        self.assertEqual(rm(src), 0)

    def test_mkdir_rmdir(self):
        """Test directory creation, cleanup and copying.
        """
        dest = os.path.join(tempfile.gettempdir(), 'ixpeobssim_temp')
        if os.path.exists(dest):
            self.assertEqual(rmdir(dest), 0)
        self.assertEqual(mkdir(dest), 0)
        self.assertEqual(rmdir(dest), 0)


if __name__ == '__main__':
    unittest.main()

#!/urs/bin/env python
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

"""Interface to the HXRB catalog described in Liu et al., 2006
https://www.aanda.org/articles/aa/abs/2006/33/aa4987-06/aa4987-06.html
"""

from __future__ import print_function, division

import os

from astropy.io import fits

from ixpeobssim import IXPEOBSSIM_TARGETS_DATA


class xHMXBCatalog:

    """
    """

    def __init__(self):
        """
        """
        file_path = os.path.join(IXPEOBSSIM_TARGETS_DATA, 'J_A+A_455_1165_table1.dat.fits')
        self._table = fits.open(file_path)[1]
        print(self._table)



if __name__ == '__main__':
    cat = xLMXBCatalog()

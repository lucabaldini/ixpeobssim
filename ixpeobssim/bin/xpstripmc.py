#!/usr/bin/env python
#
# Copyright (C) 2021 the ixpeobssim team.
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


__description__ = \
"""Remove all the Monte Carlo information from ixpeobssim photon lists.
"""

from astropy.io import fits

from ixpeobssim.core.fitsio import read_hdu_list_in_memory
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger


PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_suffix(default='nomc')
PARSER.add_overwrite()

_MC_EXTENSIONS = ('MONTE_CARLO', 'ROITABLE')

def xpstripmc(**kwargs):
    """
    """
    outlist = []
    for file_path in kwargs.get('filelist'):
        hdu_list = read_hdu_list_in_memory(file_path)
        hdu_list = fits.HDUList([hdu for hdu in hdu_list if hdu.name not in _MC_EXTENSIONS])
        outfile = file_path.replace('.fits', '_%s.fits' % kwargs.get('suffix'))
        logger.info('Writing output stripped HDU list to %s...', outfile)
        hdu_list.writeto(outfile, overwrite=kwargs.get('overwrite'))
        outlist.append(outfile)
    return outlist



def main():
    xpstripmc(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

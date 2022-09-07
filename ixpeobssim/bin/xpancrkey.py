#!/usr/bin/env python
#
# Copyright (C) 2022, the ixpeobssim team.
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

"""xpancrkey app.
"""

from __future__ import print_function, division

from astropy.io import fits

from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger, abort


#pylint: disable=invalid-name, too-many-locals, no-member


__description__ = \
"""Update the ANCRFILE keyword in the primary header of the target binned file.

This is germane to what the fparkey FTOOL
https://heasarc.gsfc.nasa.gov/lheasoft/ftools/fhelp/fparkey.html
is doing, except for the much more limited scope. On the plus side, this small
app will perform some basic check on the value being written in the header in
order to avoid common pitfalls.
"""

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_argument('--arffiles', type=str, nargs='+', default=[],
    help='path(s) to the new .arf files to use.')
PARSER.add_suffix('ancrkey')
PARSER.add_overwrite()


EXTENSION_DICT = {'PHA1': 'arf', 'PHA1Q': 'mrf', 'PHA1U': 'mrf'}
SUPPORTED_BIN_ALGS = tuple(EXTENSION_DICT.keys())


def update_header(file_path, arf_file_path, **kwargs):
    """Update the ANCRFILE keyword in a binned file.
    """
    hdu_list = fits.open(file_path)
    bin_alg = hdu_list['PRIMARY'].header['BINALG']
    if bin_alg not in EXTENSION_DICT:
        abort('Only %s files are supported, not %s' % (SUPPORTED_BIN_ALGS, bin_alg))
    if not arf_file_path.endswith(EXTENSION_DICT[bin_alg]):
        abort('Target ANCRFILE %s has the wrong extension (%s expected)' % \
        (arf_file_path, EXTENSION_DICT[bin_alg]))
    logger.info('Updating ANCRFILE keyword of %s to %s...', file_path, arf_file_path)
    hdu_list['SPECTRUM'].header['ANCRFILE'] = arf_file_path
    outfile = file_path.replace('.fits', '_%s.fits' % kwargs.get('suffix'))
    logger.info('Writing modified file to %s...', outfile)
    hdu_list.writeto(outfile, overwrite=kwargs.get('overwrite'))
    return outfile


def xpancrkey(**kwargs):
    """Run the app.
    """
    file_list = kwargs.get('filelist')
    arf_files = kwargs.get('arffiles')
    if len(file_list) != len(arf_files):
        abort('The lengths of the file lists do not match')
    return [update_header(*args, **kwargs) for args in zip(file_list, arf_files)]


def main():
    """main() entry point.
    """
    xpancrkey(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()

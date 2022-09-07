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

import os

import ixpeobssim.utils.system_ as system_
from ixpeobssim.utils.logging_ import logger


__description__ = \
"""Small Python wrapper around the HEASARC GRPPHA utility, see
https://heasarc.gsfc.nasa.gov/ftools/caldb/help/grppha.txt

In a nutshell, you should be able to pass any command that you would pass to
grppha interactively via the --comm command-line switch, e.g.

> xpgrppha.py --comm "group min 100" pha1.fits

will (loosely) map to:

> grppha infile="pha1.fits" outfile="pha1_grppha.fits" comm="group min 100 & write & exit" chatter=5 clobber=yes


(Note that "& write & exit" are automatically added at the end.)

While this is supposed to be as faithful as possible to the original, underlying
application, there are some notable differences you should keep in mind:

* you can provide an arbitrary number of input files as command-line
  arguments, and the wrapper will happily iterate over them;
* sticking to the ixpeobssim convensions, the name of the output file is
  programmatically built from the input, and you can (at least partially)
  control that via the --suffix command-line switch;
* the wrapper, as all the other ixpeobssim applications, offers a --overwrite
  command-line option that is similar in spirit to the grppha clobber
  option (for technical reasons grppha is always run with the clobber option
  set to "yes" and the check on the output file is done independently).
"""


def xpgrppha(**kwargs):
    """Wrapper implementation.
    """
    file_list = kwargs.get('filelist')
    suffix = kwargs.get('suffix')
    outlist = []
    for in_file_path in file_list:
        assert in_file_path.endswith('.fits')
        out_file_path = in_file_path.replace('.fits', '_%s.fits' % suffix)
        if not kwargs.get('overwrite') and os.path.exists(out_file_path):
            logger.warning('Output file %s already exists, skipping...', out_file_path)
            logger.warning('Use the "--overwrite True" option to overwrite')
            continue
        args = in_file_path, out_file_path, kwargs.get('comm'), kwargs.get('chatter')
        cmd = 'grppha infile="%s" outfile="%s" '\
              'comm="%s & write & exit" chatter=%d clobber=yes' % args
        system_.cmd(cmd, verbose=True)
        outlist.append(out_file_path)
    return outlist



"""Command-line switches.
"""
from ixpeobssim.utils.argparse_ import xArgumentParser

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_suffix('grppha')
PARSER.add_argument('--comm', type=str, required=True,
                    help='the GRPPHA command string')
PARSER.add_argument('--chatter', type=int, default=5,
                    help='value of the GRPPHA chatter commad-line switch')
PARSER.add_overwrite()



def main():
    args = PARSER.parse_args()
    xpgrppha(**args.__dict__)



if __name__ == '__main__':
    main()

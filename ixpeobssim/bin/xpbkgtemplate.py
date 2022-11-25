#!/usr/bin/env python
#
# Copyright (C) 2016--2022, the ixpeobssim team.
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
"""
Create a background template model starting from a series of PHA1 background files.

The PHA1 files should be prepared from a dark field with a arbitrary shapes that
have been cut from one or more files and have a BACKSCAL keyword defined.

This is suitable both for residual and total background, depending on the user
needs.

The count spectrum, normalized by the backscal, the fiducial area of the
detector and by the bin width, is parametrized with a non interpolated
spline and written to file to be used later with an appropriate config file.

The output file is written on a regular energy grid as a simple text file
with two columns---energy and background rate.
"""

import os


from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.bkg import instr
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger


PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_argument('--ssmooth', type=float, default=5.e-5,
        help='The smoothing coefficient ("s" argument in the scipy documentation) \
        used for the non interpolating spline. Note this is very important, as \
        it controls the level at which the spline is capturing the fluctuations \
        of the input data points. (s=0 is effectively an interpolating spline, \
        but the actual value depends on the scale of the input data and it is \
        not trivial to establish a priori.)')
PARSER.add_outfile(default=os.path.join(IXPEOBSSIM_SRCMODEL, 'ascii',
        'instrumental_bkg_template.txt'))


def xpbkgtemplate(**kwargs):
    filelist = kwargs.get('filelist')
    ssmooth = kwargs.get('ssmooth')
    outfile = kwargs.get('outfile')
    logger.info (f'loading background spectra from {filelist}...')
    instr.create_backgound_template(filelist, ssmooth=ssmooth, outfile=outfile)



def main():
    xpbkgtemplate(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

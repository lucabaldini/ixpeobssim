#!/usr/bin/env python
#
# Copyright (C) 2016--2020, the ixpeobssim team.
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
"""Generic display interface to ixpeobssim binned data products.

This application allows to display the content of xpbin output files in all
the possible flavors. It is equipped to recognize automagically the binning
algotithm used to create the FITS file and use the proper facilty for the
display, i.e., you should be able to pass any binned product coming out of
xpbin and xpbinview should do the right thing to display that file to you.
"""

import os

from astropy.io import fits

from ixpeobssim.binning import read_binned_file_list
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.utils.logging_ import logger, abort


def xpbinview(**kwargs):
    """Quick binned-file viewer interface.
    """
    file_list = kwargs.get('filelist')
    outfile = kwargs.get('outfile')
    try:
        bin_alg = fits.open(file_list[0])[0].header['BINALG']
    except Exception as e:
        abort('Could not determine file type (%s)' % e)
    viewer = read_binned_file_list(bin_alg, file_list)
    logger.info(viewer)

    _kwargs = dict()
    if bin_alg in ('PHA1', 'PHA1Q', 'PHA1U'):
        _kwargs['ylabel'] = '%s counts [s$^{-1}$]' % bin_alg
    elif bin_alg in ('PHA1QN', 'PHA1UN'):
        _kwargs['ylabel'] = '%s counts' % bin_alg
    elif bin_alg in ('CMAP',):
        _kwargs['stretch'] = kwargs['stretch']
        _kwargs['zlabel'] = kwargs.get('zlabel') or 'Counts/pixel'
        vmin, vmax = kwargs.get('vmin'), kwargs.get('vmax')
        if vmin is not None:
            _kwargs['vmin'] = vmin
        if vmax is not None:
            _kwargs['vmax'] = vmax
    elif bin_alg in ('LC',):
        _kwargs['mjd'] = kwargs['mjd']
    viewer.plot(**_kwargs)
    if outfile is not None:
        if os.path.exists(outfile) and not kwargs.get('overwrite'):
            logger.info('Output file %s already exists.', outfile)
            logger.info('Remove it or set "overwrite = True" to overwite it.')
        else:
            logger.info('Saving current figure to %s...', outfile)
            plt.savefig(outfile, dpi=kwargs.get('dpi'))
    plt.show()


# Command-line switches.
PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_outfile()
PARSER.add_overwrite()
PARSER.add_stretch()
PARSER.add_vrange()
PARSER.add_argument('--zlabel', type=str, default=None,
    help='the color map laber for count maps')
PARSER.add_argument('--dpi', type=int, default=100,
    help='resolution of the output image in dot per inches')
PARSER.add_boolean('--mjd', default=False,
                  help='Show light curves in MJD time (as opposed to MET)')


def main():
    xpbinview(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()

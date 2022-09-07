#!/usr/bin/env python
#
# Copyright (C) 2016, 2018, the ixpeobssim team.
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

"""xpselect app.
"""

from __future__ import print_function, division

import os

from ixpeobssim.utils.logging_ import logger
from ixpeobssim.evt.subselect import xEventSelect
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.os_ import check_input_file


__description__ = \
"""Basic data selection interface. This utility can process IXPE event files
(i.e., photon lists) and apply a generic selection on sky-position, time,
phase, energy, etc., producing smaller event files. In conjunction with
xpbin this is handy to, e.g., energy- or phase-resolved polarization
measurements.

xpselect supports spatial selections from a ds9 region file thorough the
--regfile command-line switch.
"""

LIVETIME_UPDATE_ALGS = ('LTSUM', 'LTSCALE')
LIVETIME_UPDATE_DEFAULT = 'LTSCALE'

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_suffix()
PARSER.add_tbounds(None, None)
PARSER.add_boolean('--tinvert', default=False,
                   help='invert the time selection')
PARSER.add_phasebounds(None, None)
PARSER.add_boolean('--phaseinvert', default=False,
                   help='invert the phase selection')
PARSER.add_ebounds(None, None)
PARSER.add_boolean('--einvert', default=False,
                   help='invert the energy selection')
PARSER.add_argument('--ra', type=float, default=None,
                    help='RA of acceptance cone in decimal degrees')
PARSER.add_argument('--dec', type=float, default=None,
                    help='Dec of acceptance cone in decimal degrees')
PARSER.add_argument('--rad', type=float, default=None,
                    help='ROI radius in arcminutes')
PARSER.add_argument('--innerrad', type=float, default=None,
                    help='ROI inner radius in arcminutes (for selecting annuli)')
PARSER.add_argument('--regfile', type=str, default=None,
                    help='path to a ds9 region file')
PARSER.add_boolean('--reginvert', default=False,
                   help='invert the ds9 region file selection')
PARSER.add_argument('--mask', type=str, default=None,
                    help='path to .npy file with a saved boolean event mask')
PARSER.add_argument('--mcsrcid', action='append', type=int, default=[],
                    help='the Monte Carlo source ID to select')
PARSER.add_mc()
PARSER.add_boolean('--ltimeupdate', default=False,
                   help='update the livetime-related keywords in the output files')
PARSER.add_argument('--ltimealg', choices=LIVETIME_UPDATE_ALGS, default=LIVETIME_UPDATE_DEFAULT,
                   help='algorithm to be used to update the livetime')
PARSER.add_overwrite()



def xpselect(**kwargs):
    """Application for data subselection.

    We want to (loosely) model this on
    http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtselect.txt
    """
    
    file_list = kwargs.get('filelist')
    outlist = []
    for file_path in file_list:
        check_input_file(file_path, 'fits')
        event_select = xEventSelect(file_path, **kwargs)
        outfile = event_select.get('outfile')
        outlist.append(outfile)
        if os.path.exists(outfile) and not event_select.get('overwrite'):
            logger.info('Output file %s already exists.', outfile)
            logger.info('Remove it or set "overwrite = True" to overwite it.')
        else:
            event_select.select()
    return outlist


def main():
    """main() entry point.
    """
    xpselect(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

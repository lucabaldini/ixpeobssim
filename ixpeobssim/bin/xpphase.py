#!/usr/bin/env python
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

from __future__ import print_function, division


__description__ = \
"""Utility to build phase column into event file.

The program generates a phase array using pulsar time and ephemeris
and builds a new event FITS file with the new PHASE column.
"""

import os

from astropy.io import fits

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.srcmodel.ephemeris import xEphemeris


def _get_ephemeris(**kwargs):
    """Select the chosen ephemeris setup.

    We pass the ephemeris by either a par file or specifying the actual
    numbers one by one. Since argparse does not support out of the box
    mutually exclusive complex groups, see, e.g.
    https://stackoverflow.com/questions/17909294
    we pefer to encapsulate the small bit of logic involved in sorting
    out which method we are using in this function, rather than going for
    subcommands.

    That all said, the logic is: if a par file is passed we use that and
    make sure that none of the other mutuallt exclusive options are passed.
    Otherwise we go ahead and build an ephemeris from the numbers.

    Warning
    -------
    I don't really like using the PARSER global variable in here, but emitting a
    parser error when the arguments are inconsistent seems like the most reasonable
    thing to do, and I don't think it's worth subclassing the argument
    parser and adding a specialize parse_arguments() slot for this.
    """
    ephem_keys = ['met0', 'nu0', 'nudot0', 'nuddot']
    par_file = kwargs.get('parfile')
    if par_file is not None:
        for key in ephem_keys:
            if kwargs.get(key) is not None:
                PARSER.error('Cannot specify --%s when using the --parfile option' % key)
        if not par_file.endswith('.par'):
            PARSER.error('The --parfile switch should point to a .par file')
        return xEphemeris.from_file(par_file)
    # We need at the very least n0 to do something sensible.
    if kwargs.get('nu0') is None:
        PARSER.error('At least --nu0 should be specified')
    # And now we set to 0. the parameters that are still None. (Remember, we
    # use None to signal the the command-line switch have not been set, but
    # for these three parameters we have a sensible default value of 0.)
    for key in ['met0', 'nudot0', 'nuddot']:
        if kwargs[key] is None:
            kwargs[key] = 0.
    params = {key: kwargs.get(key) for key in ephem_keys}
    return xEphemeris(**params)


def xpphase(**kwargs):
    """Build the phase column into event file.
    """
    file_list = kwargs.get('filelist')
    eph = _get_ephemeris(**kwargs)

    outlist = []
    for file_path in file_list:
        assert file_path.endswith('.fits')
        if not os.path.exists(file_path):
            logger.warning('Input file %s does not exist!', file_path)
            logger.warning('Check your input file.')
            continue

        suffix = kwargs.get('suffix')
        #Path to FITS file with phase
        outfile = file_path.replace('.fits', '_%s.fits' % suffix)
        outlist.append(outfile)
        if os.path.exists(outfile) and not kwargs.get('overwrite'):
            logger.info('Output file %s already exists.', outfile)
            logger.info('Remove the file or set "overwrite = True" to overwrite it.')
            continue

        logger.info('Opening "%s"...', file_path)
        hdu = fits.open(file_path)
        evt_hdu = hdu['EVENTS']
        time_ = evt_hdu.data['TIME']
        evt_header = evt_hdu.header

        logger.info('Calculating pulsar phase...')
        phase = eph.fold(time_, evt_header['TSTART'], kwargs.get('phi0'))
        logger.info('Creating phase column...')
        col_ = fits.Column(name='PHASE', array=phase, format='E')
        new_hdu = fits.BinTableHDU.from_columns(evt_hdu.data.columns + col_, header=evt_header)

        # Fix missed comments in new header
        ncols = len(evt_hdu.data.columns)
        for i in range(ncols):
            comment = evt_header.comments['TTYPE%d' % (i + 1)]
            new_hdu.header.comments['TTYPE%d' % (i + 1)] = comment
        new_hdu.header.comments['TTYPE%d' % (ncols + 1)] = 'event phase (periodic sources)'

        logger.info('Writing to %s...', outfile)
        hdulist = fits.HDUList([hdu[0], new_hdu] + hdu[2:])
        hdulist.writeto(outfile, overwrite=True)
        hdulist.info()
        hdu.close()
        hdulist.close()
    logger.info('Done!')
    return outlist



"""Command-line switches.
"""
from ixpeobssim.utils.argparse_ import xArgumentParser

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_suffix('phase')
PARSER.add_argument('--parfile', type=str, default=None,
                    help='path to the input parameter file')
PARSER.add_argument('--met0', type=float, default=None,
                    help='reference MET of the ephemeris in s')
PARSER.add_argument('--nu0', type=float, default=None,
                    help='frequency at t0 in Hz')
PARSER.add_argument('--nudot0', type=float, default=None,
                    help='time first derivative of frequency at t0 in 1/(s^2)')
PARSER.add_argument('--nuddot', type=float, default=None,
                    help='time second derivative of frequency in 1/(s^3)')
PARSER.add_phi0()
PARSER.add_overwrite()



def main():
    xpphase(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

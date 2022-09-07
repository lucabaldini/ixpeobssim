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
"""Utility to build orbital phase column into event file.
"""

import os

import numpy
from astropy.io import fits

from ixpeobssim import  IXPEOBSSIM_DATA
from ixpeobssim.srcmodel.ephemeris import xOrbitalEphemeris
from ixpeobssim.srcmodel.tdelays import xTDelays
from ixpeobssim.utils.logging_ import logger


def barycorr(T, ephemeris, **kwargs):
    """Apply all time correction available.
    """
    _times = xTDelays(T, unit='met', name='TIME')
    _times.apply_delay(ephemeris, ra=kwargs['ra'], dec=kwargs['dec'], delay='all', binary=True)
    return _times.metvalue

def get_ophase(T, ephemeris, phi0=0.):
    """Calculate the event orbital phase given event time and source parameter file.
    """
    dt = (T - ephemeris.t_orbital)
    phase = phi0 + dt * ephemeris.omega_orb(dt) # PBDOT in
    ph = phase - numpy.floor(phase)
    return numpy.array(ph)

def xpophase(**kwargs):
    """Build the phase column into event file and (if required) update the TIME
    column with barycentred times.
    """
    file_list = kwargs.get('filelist')

    ephemeris = xOrbitalEphemeris.from_file(kwargs['parfile'])

    outlist = []
    for file_path in file_list:
        assert file_path.endswith('.fits')
        if not os.path.exists(file_path):
            logger.warning('Input file %s does not exist!', file_path)
            logger.warning('Check your input file.')
            continue

        suffix = kwargs.get('suffix')
        # Path to FITS file with phase
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

        if kwargs.get('bary'):
            # Update TIME data column and header comment
            logger.info('Applying barycentric time correction...')
            time_ = barycorr(time_, ephemeris, **kwargs)
            evt_hdu.data['TIME'] = time_
            evt_header.comments['TTYPE2'] = 'event barycentred time in seconds'

        logger.info('Calculating orbital phase...')
        ophase = get_ophase(time_, ephemeris)
        logger.info('Creating phase column...')
        col_ = fits.Column(name='ORB_PHASE', array=ophase, format='E')
        new_hdu = fits.BinTableHDU.from_columns(evt_hdu.data.columns + col_,
                                                header=evt_header)

        # Fix missed comments in new header
        ncols = len(evt_hdu.data.columns)
        for i in range(ncols):
            comment = evt_header.comments['TTYPE%d' % (i + 1)]
            new_hdu.header.comments['TTYPE%d' % (i + 1)] = comment
        new_hdu.header.comments['TTYPE%d' % (ncols + 1)] = 'event orbital phase (binary sources)'

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
import ast
from ixpeobssim.utils.argparse_ import xArgumentParser

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_suffix('ophase')
PARSER.add_argument('--parfile', type=str, required=True,
                    help='path to the input parameter file')
PARSER.add_argument('--ra', type=float, default=None,
                    help='source right ascension in deg')
PARSER.add_argument('--dec', type=float, default=None,
                    help='source declination in deg')
PARSER.add_phi0()
PARSER.add_argument('--bary', type=ast.literal_eval,
                    default=False, choices=[True, False],
                    help='apply or not barycentric time correction')
PARSER.add_overwrite()


def main():
    xpophase(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

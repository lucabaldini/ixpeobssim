#!/usr/bin/env python
#
# Copyright (C) 2015--2020, the ixpeobssim team.
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
"""Basic event binning interface.

This application allows to bin event files (i.e., photon lists) in a fairly
generic fashion, creating count spectra, light curves, maps and more advances
data products such as those for storing the Stokes parameters.
"""


import os
import sys

from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.os_ import check_input_file
from ixpeobssim.binning import BINNING_WRITE_DICT


BIN_ALGS = list(BINNING_WRITE_DICT.keys())
TBIN_ALGS = ['FILE', 'LIN', 'LOG']
EBIN_ALGS = ['FILE', 'LIN', 'LOG', 'EQP', 'LIST']
PRJCTS = ['AIT', 'ZEA', 'ARC', 'CAR', 'GLS', 'MER', 'NCP', 'SIN', 'STG', 'TAN']
COORD_SYS = ['CEL', 'GAL']
UNIVERSAL_KWARGS = ['algorithm', 'help', 'overwrite', 'suffix']

# Complete the help based on the actual implementation of the binning classes.
__description__ += '\nSupported binning algorithms\n'
for alg, class_ in BINNING_WRITE_DICT.items():
    opts = ', '.join(class_.SUPPORTED_KWARGS + ['overwrite', 'suffix'])
    __description__ += '* %s: %s\n   --%s\n' % (alg, class_.INTENT, opts)


def _check_kwargs(**kwargs):
    """Small facility to enforce that the keyword arguments passed from
    command-line make sense for the specific binning algorithm.

    This implies that the underlying event binning classes exposes a
    top-level SUPPORTED_KWARGS overriding the default of the xEventBinningBase
    class (None, which means the check is disengaged).

    See https://bitbucket.org/ixpesw/ixpeobssim/issues/328
    """
    algorithm = kwargs['algorithm']
    binning_class = BINNING_WRITE_DICT[algorithm]
    supported_keys = binning_class.SUPPORTED_KWARGS
    # If the binning class does not provide a list of supported keys, we
    # proceed silently.
    if supported_keys is None:
        return
    # Add the command-line switches that are supported by all algoritms to the
    # list of supported keys.
    supported_keys += UNIVERSAL_KWARGS
    # Parse the command line switches from sys.argv
    cmd_line_switches = [item.strip('-') for item in sys.argv if item.startswith('--')]

    def _find_key(switch_, kwargs_keys):
        """Small convenience function to match partial switch names to actual
        keys in the argument parser namespace.

        This takes care of the fact that, e.g., "--algorithm" can be spelled
        as "--alg" from command-line.
        """
        if switch_ in kwargs_keys:
            return switch_
        for key in kwargs_keys:
            if key.startswith(switch_):
                return key
        logger.error('Could not map "%s" to a known xpbin option.', switch_)
        logger.error('sys.argv: %s', sys.argv)
        return switch_

    # Make sure we have all the command line switches in the complete, long form.
    kwargs_keys = list(kwargs.keys())
    cmd_line_switches = [_find_key(item, kwargs_keys) for item in cmd_line_switches]
    logger.info('Checking command-line options: %s...', cmd_line_switches)
    for key in cmd_line_switches:
        if key not in supported_keys:
            # Oopss... we have an unsupported command-line switch and we do exit.
            PARSER.print_help()
            logger.error('Option "%s" not supported for "%s" binning algorithm.', key, algorithm)
            abort()


def xpbin(**kwargs):
    """Application to bin the data.

    We want to (loosely) model this on
    http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtbin.txt
    """
    #_check_kwargs(**kwargs)
    outlist = []
    for file_path in kwargs.get('filelist'):
        check_input_file(file_path, 'fits')
        event_binning = BINNING_WRITE_DICT[kwargs['algorithm']](file_path, **kwargs)
        outfile = event_binning.get('outfile')
        outlist.append(outfile)
        if os.path.exists(outfile) and not event_binning.get('overwrite'):
            logger.info('Output file %s already exists.' % outfile)
            logger.info('Remove it or set "overwrite = True" to overwite it.')
        else:
            event_binning.bin_()
    return outlist


"""Command-line switches.
"""
from ixpeobssim.utils.argparse_ import xArgumentParser

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_argument('--algorithm', choices=BIN_ALGS, required=True,
                    help='the binning algorithm')
PARSER.add_suffix(None)
PARSER.add_irfname(default=None)
PARSER.add_grayfilter()
PARSER.add_boolean('--acceptcorr', default=True,
    help='enable/disable the acceptance correction for polarization cubse and maps')
PARSER.add_weights(default=False)
PARSER.add_weightcol()
PARSER.add_argument('--tbinalg', choices=TBIN_ALGS, default='LIN',
                    help='time binning specification')
PARSER.add_tbounds(None, None)
PARSER.add_argument('--tbins', type=int, default=100,
                    help='number of bins for LIN/LOG time binning')
PARSER.add_argument('--tbinfile', type=str, default=None,
                    help='path to the optional time bin definition file')
PARSER.add_argument('--phasebins', type=int, default=50,
                    help='number of bins for phase binning')
PARSER.add_argument('--npix', type=int, default=200,
                    help='number of pixels per side in the output image')
PARSER.add_argument('--pixsize', type=float, default=None,
                    help='the pixel size of the output image in arcseconds')
PARSER.add_argument('--xref', type=float, default=None,
                    help='the horizontal position of the image center')
PARSER.add_argument('--yref', type=float, default=None,
                    help='the vertical position of the image center')
PARSER.add_argument('--proj', choices=PRJCTS, default='TAN',
                    help='coordinate projection')
PARSER.add_ebinning(default_emin=2., default_emax=8., ebin_algs=EBIN_ALGS,
                    default_ebinalg='LOG', default_ebins=4)
PARSER.add_argument('--phibins', type=int, default=75,
                    help='number of bins for LIN/LOG phi binning')
PARSER.add_argument('--thetabins', type=int, default=17,
                    help='number of bins for off-axis angle binning')
PARSER.add_mc()
PARSER.add_overwrite()


def main():
    xpbin(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

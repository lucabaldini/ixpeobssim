#!/usr/bin/env python
#
# Copyright (C) 2018--2020, the ixpeobssim team.
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
"""Minimal interface to perform a spectro-polarimetric fit in XSPEC, via the
Python bindings.

This is primarily designed to fit a spectral and polarimetric model to a series
of I, Q and U Stokes spectra (typically nine input files---three spectra for
each detector unit) or a purely polarimetric fit to a series of normalized
Q/I and U/I spectra (typically six input files).

Although only a small part of the flexibility provided by a real XSPEC session
is currently exposed, this simple application is a sensible mean to get started
with the use of the ixpeobssim polarimetric local models within XSPEC.

Note that the application fully supports purely spectral or polarimetric fits
in a transparent fashion.
"""

import os

from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.utils.logging_ import logger, abort

from ixpeobssim.utils.environment import PYXSPEC_INSTALLED
if PYXSPEC_INSTALLED:
    import ixpeobssim.evt.xspec_ as xspec_
    xspec_stat_methods = xspec_.STAT_METHODS
    xspec_.load_local_models()
else:
    xspec_stat_methods = []



def xpxspec(**kwargs):
    """Do a spectro-polarimetric fit in XSPEC
    """
    if not PYXSPEC_INSTALLED:
        return
    file_list = kwargs.get('filelist')
    num_files = len(file_list)
    xspec_.load_input_files(*file_list)
    xspec_.select_energy_range(kwargs.get('emin'), kwargs.get('emax'))
    xspec_.setup_fit_model(kwargs.get('model'), kwargs.get('params'))
    for par_index, par_value in kwargs.get('fixpars'):
        logger.info('Setting parameter %d to %f...', par_index, par_value)
        xspec_.fix_parameter(par_index, par_value)
    if kwargs['fit']:
        args = kwargs.get('statmethod'), kwargs.get('niterations'), kwargs.get('error')
        fit_output = xspec_.fit(*args)
        logger.info(fit_output)
    if kwargs['plot']:
        xspec_.plot(kwargs.get('figname'))
        plt.show()
    return fit_output




"""Command-line switches.
"""
import ast

def fixpars_opt(arg):
    """Custom formatter for the fixpars command-line switch.
    """
    par_id, par_val = arg.split(',')
    return int(par_id), float(par_val)


from ixpeobssim.utils.argparse_ import xArgumentParser

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_argument('--model', type=str, required=True,
                    help='the spectral model for the fit')
PARSER.add_argument('--params', type=ast.literal_eval, default=None,
                    help='the initial values of the fit parameters')
PARSER.add_argument('--fixpars', type=fixpars_opt, nargs='+', default=[],
                    help='freeze one or more fit-parameter value(s)')
PARSER.add_argument('--statmethod', type=str, default='chi',
                    choices=xspec_stat_methods,
                    help='the fit statistics to be used')
PARSER.add_argument('--niterations', type=int, default=25,
                    help='maximum  number of fit iterations')
PARSER.add_ebounds()
PARSER.add_boolean('--fit', default=True,
                   help='fit the data')
PARSER.add_boolean('--error', default=True,
                   help='calculate confidence intervals for fit parameters')
PARSER.add_boolean('--plot', default=True,
                   help='plot the fit results')
PARSER.add_argument('--figname', type=str, default='XSPEC fit',
                    help='base name for the output figures')
PARSER.add_overwrite()


def main():
    args = PARSER.parse_args()
    xpxspec(**args.__dict__)



if __name__ == '__main__':
    main()

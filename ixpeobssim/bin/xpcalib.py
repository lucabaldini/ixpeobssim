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

import os

import numpy

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.evt.gti import xSimpleGTIList
from ixpeobssim.instrument import DU_IDS, du_suffix
from ixpeobssim.irf import load_irf_set
from ixpeobssim.srcmodel import import_roi
from ixpeobssim.srcmodel.calibsrc import xCalibrationROIModel
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.profile import xChrono
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.time_ import string_to_met_utc


__description__ = \
"""Create photon lists from IXPE calibration runs, either on ground or in orbit.
"""


def xpcalib(**kwargs):
    """Run the calibration data-taking run.

    This has a lot in common with xpobssim.py, but since the latter is already
    kind of complex, accepting a bit of code duplication rathen than making
    xpobssim even more comples seemed like the right thing to do.
    """
    assert kwargs.get('configfile').endswith('.py')
    if kwargs.get('outfile') is None:
        outbase = os.path.basename(kwargs['configfile']).replace('.py', '')
        kwargs['outfile'] = os.path.join(IXPEOBSSIM_DATA, outbase)
        logger.info('Setting output file base name to %s...', kwargs['outfile'])
    outfile = kwargs['outfile'].replace('.fits', '')
    # Handle the random seed. According to the numpy documentation, the
    # random seed used to initialize the pseudo-random number generator can be
    # any integer between 0 and 2**32 - 1 inclusive (among other different
    # things that numpy is able to convert to an integer internally).
    # If the seed is None it is chosen randomly, see
    # https://bitbucket.org/ixpesw/ixpeobssim/issues/182
    random_seed = kwargs['seed']
    if random_seed is None:
        random_seed = numpy.random.randint(0, 2**32)
    chrono = xChrono()
    # Setup the source model.
    logger.info('Setting up the source model...')
    roi_model = import_roi(kwargs.get('configfile'))
    if not isinstance(roi_model, xCalibrationROIModel):
        logger.error('The source model does not seem to be a calibration setup, please run xpobssim.py')
        abort()
    logger.info(roi_model)
    logger.info('Done %s.', chrono)
    # Convert the start_date and duration into start_met and stop_met, which
    # are frequently consumed downstream.
    kwargs['start_met'] = string_to_met_utc(kwargs.get('startdate'), lazy=True)
    kwargs['stop_met'] = kwargs.get('start_met') + kwargs.get('duration')
    # Calculate the relevant GTIs.
    kwargs['gti_list'] = xSimpleGTIList(kwargs['start_met'], kwargs['stop_met'])
    outlist = []
    for du_id in DU_IDS:
        kwargs['outfile'] = '%s_%s.fits' % (outfile, du_suffix(du_id))
        outlist.append(kwargs['outfile'])
        if os.path.exists(kwargs['outfile']) and not kwargs['overwrite']:
            logger.info('Output file %s already exists.', kwargs['outfile'])
            _info = 'Remove the file or set overwrite = True to overwrite it.'
            logger.info(_info)
            continue
        # Set the random seed---this has to be different for each of the DUs
        # and we just add the DU identifier (-1) to achieve that.
        _seed = random_seed + du_id - 1
        logger.info('Setting the random seed to %d...', _seed)
        numpy.random.seed(_seed)
        # Load the response functions.
        logger.info('Loading the instrument response functions...')
        irf_set = load_irf_set(kwargs.get('irfname'), du_id)
        logger.info('Done %s.', chrono)
        # Run the actual simulation.
        logger.info('Generating the photon list...')
        event_list = roi_model.rvs_event_list(irf_set, **kwargs)
        event_list.write_fits('xpcalib.py', roi_model, irf_set, **kwargs)
        logger.info('Done for detector unit # %i %s.', du_id, chrono)
    logger.info('All done %s!', chrono)
    return outlist


"""Command-line switches."""

PARSER = xArgumentParser(description=__description__)
PARSER.add_outfile()
PARSER.add_configfile(required=True)
PARSER.add_irfname()
PARSER.add_duration()
PARSER.add_ebounds(default_emin=1., default_emax=12.)
PARSER.add_startdate()
PARSER.add_seed()
PARSER.add_charging()
PARSER.add_deadtime()
PARSER.add_boolean('--lv1a', False, help='create a (pseudo) Lv1a file')
PARSER.add_overwrite()


def main():
    xpcalib(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

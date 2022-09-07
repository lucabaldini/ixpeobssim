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

"""xppicorr app.
"""

from __future__ import print_function, division

from astropy.io import fits
import numpy

from ixpeobssim.evt.event import xEventFile, xEventList
from ixpeobssim.evt.picorr import xPulseInvariantCorrection
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.misc import process_file_list
from ixpeobssim.utils.os_ import check_input_file, check_output_file

#pylint: disable=no-member

__description__ = \
"""Correct the PI column in a level-2 file, given either a constant scale factor
and/or offset, or a proper FITS file with a time-dependent correction.

In order to avoid steps in the output corrected files, a floating-point smearing
factor in the [-0.5, 0.5] interval is added to the original PI before the
correction. When run in deterministic mode this factor is determined by the
event time (more precisely, using the fractional part of the timestamps), otherwise
we just throw a random number with a uniform distribution.

Note this is mainly intended for debugging purposes---you are changing the
actual data, so you'd better know what you are doing :-)
"""

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_argument('--offset', type=float, default=0.,
    help='global offset for the PI column')
PARSER.add_argument('--slope', type=float, default=1.,
    help='global scale factor for the PI column')
PARSER.add_argument('--corrfile', type=str, default=None,
    help='path to the FITS file for a time-dependent correction')
PARSER.add_boolean('--deterministic', default=True,
    help='run in deterministic mode')
PARSER.add_seed()
PARSER.add_suffix('picorr')
PARSER.add_overwrite()



def _smearing_delta_pi(event_file, **kwargs):
    """
    Calculate the smearing factor to be added to the PI prior to the correction,
    in order to avoid steps in the output file.

    Depending on whether xppicorr is run in deterministic mode or not, this is
    done using the fractional part of the timestamp, or an actual random number.
    See https://bitbucket.org/ixpesw/ixpeobssim/issues/595/ for more details.
    """
    if kwargs.get('deterministic'):
        logger.info('Applying deterministic smearing...')
        _, microseconds = xEventList.split_event_time(event_file.time_data())
        delta = 1.e-6 * microseconds - 0.5
    else:
        seed = kwargs.get('seed')
        if seed is not None:
            logger.info('Applying random smearing...')
            logger.info('Setting the random seed to %d...', kwargs.get('seed'))
            numpy.random.seed(seed)
        delta = numpy.random.uniform(-0.5, 0.5, size=event_file.num_events())
    logger.debug('Smearing stat (min, max, average): %.6f, %.6f, %.6f',
        delta.min(), delta.max(), delta.mean())
    return delta


def _process_file(file_path, **kwargs):
    """Process a single file.
    """
    output_file_path = check_output_file(file_path, kwargs.get('suffix'), kwargs.get('overwrite'))
    if output_file_path is None:
        return None
    event_file = xEventFile(file_path)
    # Retrieve the PI.
    pi = event_file.pi_data().astype(float)
    # A random fluctuation not to have steps in the output binned files.
    pi += _smearing_delta_pi(event_file, **kwargs)
    if kwargs.get('corrfile') is not None:
        correction = xPulseInvariantCorrection(kwargs.get('corrfile'))
        slope, offset = correction(event_file.time_data())
        logger.info('Average slope: %.3f', slope.mean())
        logger.info('Average offset: %.3f', offset.mean())
    else:
        slope, offset = kwargs.get('slope'), kwargs.get('offset')
    # Apply the correction.
    pi = pi * slope + offset
    # Convert back to integer---note we are clipping the PI values below zero, if
    # any after the correction.
    pi = numpy.rint(pi).clip(0)
    event_file.set_column('EVENTS', 'PI', pi)
    event_file.write(output_file_path, overwrite=kwargs.get('overwrite'))
    return output_file_path


def xppicorr(**kwargs):
    """Application entry point.
    """
    return process_file_list(_process_file, kwargs.get('filelist'), **kwargs)


def main():
    """main() entry point.
    """
    xppicorr(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()

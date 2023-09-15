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

from ixpeobssim.bin.xpobssim import _prepare_simulation
from ixpeobssim.instrument import DU_IDS, du_suffix
from ixpeobssim.irf import load_irf_set
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger



__description__ = \
"""Create a photon list to be fed into ixpesim.

This is a rough equivalent to xpobssim, except for the fact that the output FITS
file(s) represent a list of photons at the top of the GPD window, as opposed to
a list of photoelectron tracks. This allows the list to be fed into xpsim, where
the photons are propagated through the detector, and the complete information
for all the events that trigger, including the track images, is written to file.

In conjunction with ixpesim, this small application allows for a full end-to-end
workflow where arbitrarily complex source models can be simulated with the
ultimate fidelity.
"""

PARSER = xArgumentParser(description=__description__)
PARSER.add_outfile()
PARSER.add_configfile(required=True)
PARSER.add_irfname()
PARSER.add_duration()
PARSER.add_gti_settings()
PARSER.add_ebounds(default_emin=1., default_emax=12.)
PARSER.add_startdate()
PARSER.add_objname()
PARSER.add_seed()
PARSER.add_vignetting()
PARSER.add_dithering()
PARSER.add_grayfilter()
PARSER.add_roll()
PARSER.add_trajectory()
PARSER.add_sc_data()
PARSER.add_overwrite()


def xpphotonlist(**kwargs):
    """Run the application.
    """
    outfile, random_seed, roi_model, chrono, _ = _prepare_simulation(kwargs)
    outlist = []
    for du_id in DU_IDS:
        kwargs['outfile'] = '%s_%s_photon_list.fits' % (outfile, du_suffix(du_id))
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
        irf_set = load_irf_set(kwargs.get('irfname'), du_id, gray_filter=kwargs.get('grayfilter'))
        logger.info('Done %s.', chrono)
        # Run the actual simulation.
        logger.info('Generating the photon list...')
        photon_list = roi_model.rvs_photon_list(irf_set, **kwargs)
        photon_list.write_fits('xpobssim.py', roi_model, irf_set, **kwargs)
        logger.info('Done for detector unit # %d %s.', du_id, chrono)
    logger.info('All done %s!', chrono)
    return outlist


def main():
    """
    """
    xpphotonlist(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

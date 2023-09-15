#!/usr/bin/env python
#
# Copyright (C) 2015--2020 the ixpeobssim team.
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

from ixpeobssim import IXPEOBSSIM_DATA, IXPEOBSSIM_OBSDATA, IXPEOBSSIM_CALDB
from ixpeobssim.evt.event import xEventList
from ixpeobssim.instrument import DU_IDS, NUM_DETECTOR_UNITS, du_suffix
from ixpeobssim.instrument.fcw import xOnOrbitCalibrationPattern
from ixpeobssim.instrument.traj import xObservationTimeline
from ixpeobssim.irf import load_irf_set
from ixpeobssim.srcmodel import import_roi
from ixpeobssim.srcmodel.calibsrc import xCalC, xCalibrationROIModel
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.profile import xChrono
from ixpeobssim.utils.time_ import string_to_met_utc
from ixpeobssim.core.fitsio import read_hdu_list_in_memory


__description__ = \
"""Run the ixpeobssim fast simulator.

This the main application in the package, and produces a set of event files
(a.k.a. photon lists), given a source model and a set of instrument response
functions, for a given observation time.

Internally the source spectral, temporal, morphological and polarimetric
characteristics are convolved with the instrument response functions to produce
a realistic, simulated IXPE observation.
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
PARSER.add_charging()
PARSER.add_deadtime()
PARSER.add_roll()
PARSER.add_trajectory()
PARSER.add_on_orbit_calibration()
PARSER.add_timeline()
PARSER.add_sc_data()
PARSER.add_boolean('--lv1a', False, help='create a pseudo-Lv1a file')
PARSER.add_argument('--lv1version', type=int, default=5,
    help='version number for the pseudo-Lv1a support')
PARSER.add_overwrite()



def _build_timeline(roi_model, **kwargs):
    """Compile the timeline for the observation,
    """
    logger.info('Compiling the observation timeline...')
    args = kwargs.get('start_met'), kwargs.get('stop_met'), roi_model.ra,\
        roi_model.dec, kwargs.get('saa'), kwargs.get('occult')
    timeline = xObservationTimeline(*args)
    # Calculate the GTI list.
    logger.info('Assembling GTI list...')
    args = kwargs.get('gtiminduration'), kwargs.get('gtistartpad'),\
        kwargs.get('gtistoppad')
    gti_list = timeline.gti_list(*args)
    logger.info(gti_list)
    # Calculate the OCTI list.
    if kwargs.get('onorbitcalib'):
        args = kwargs.get('onorbitcalminduration'), kwargs.get('onorbitcalstartpad'),\
            kwargs.get('onorbitcalstoppad')
        octi_list = timeline.octi_list(*args)
    else:
        octi_list = []
    calib_source = xCalC(kwargs.get('onorbitcalrate'))
    args = octi_list, calib_source, kwargs.get('onorbitcaldemult')
    calib_pattern = xOnOrbitCalibrationPattern(*args)
    logger.info(calib_pattern)
    return timeline, gti_list, calib_pattern


def _map_charging_file_list(file_list, default_mapping, ext_name):
    """Tranform a charging file list (be it for charge maps or charging params)
    into a dictionary indexed by DU ID, using the default values if the file
    list is None.

    Note that we make sure that the DETNAM keyword in the proper FITS files is
    consistent with the expected DU ID, as filenames are prone to be modified
    by the user. (Since the charging map files are small, opening them twice
    should not pose a serious perfomance issue.)

    Arguments
    ---------
    file_list : 3-element list or None
        The input file list.

    default_mapping : callable
        A function mapping a DU ID into the corresponding default file path.
    """
    if file_list is None:
        file_list = [default_mapping(du_id) for du_id in DU_IDS]
    else:
        assert len(file_list) == NUM_DETECTOR_UNITS
    file_dict = dict(zip(DU_IDS, file_list))
    for du_id, file_path in file_dict.items():
        hdu_list = read_hdu_list_in_memory(file_path)
        if int(hdu_list[ext_name].header['DETNAM'][-1]) != du_id:
            abort('Inconsistent DU ID for file %s' % file_path)
    return file_dict


def _update_charging_settings(kwargs):
    """Go over the charging-related setup.

    Here we essentially turn the charging map and charging param file lists into
    dictionaries indexed by DU ID.

    Note that we have to pass the kwargs object directly, as opposed to the **variant,
    because we want to modify it in place.
    """
    # If the charging is disabled, there is nothing to do.
    if not kwargs.get('charging', False):
        return
    # Setup the charge map files.
    _folder = os.path.join(IXPEOBSSIM_OBSDATA, 'chrgmaps')
    _map = lambda du_id: os.path.join(_folder, 'ixpe_vanilla_d%d_chrgmap.fits' % du_id)
    kwargs['chrgmaps'] = _map_charging_file_list(kwargs.get('chrgmaps'), _map, 'CHRG_MAP')
    # Setup the charging param files.
    _folder = os.path.join(IXPEOBSSIM_CALDB, 'ixpe', 'gpd', 'bcf', 'chrgparams')
    _map = lambda du_id: os.path.join(_folder, 'ixpe_vanilla_d%d_chrgparams.fits' % du_id)
    kwargs['chrgparams'] = _map_charging_file_list(kwargs.get('chrgparams'), _map, 'CHRG_PAR')


def _prepare_simulation(kwargs):
    """Do all the preliminary prep-work, all the way up to the loop over the
    detector units.

    Note that we have to pass the kwargs object directly, as opposed to the **variant,
    because we want to modify it in place.
    """
    assert kwargs.get('configfile').endswith('.py')
    if kwargs.get('outfile') is None:
        outbase = os.path.basename(kwargs['configfile']).replace('.py', '')
        kwargs['outfile'] = os.path.join(IXPEOBSSIM_DATA, outbase)
        logger.info('Setting output file base name to %s...', kwargs['outfile'])
    outfile = kwargs['outfile'].replace('.fits', '')
    if kwargs['objname'] is None:
        kwargs['objname'] = os.path.basename(kwargs.get('configfile')).replace('.py', '')
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
    if isinstance(roi_model, xCalibrationROIModel):
        logger.error('The source model seems to be a calibration setup, please run xpcalib.py')
        abort()
    logger.info(roi_model)
    logger.info('Done %s.', chrono)
    # Convert the start_date and duration into start_met and stop_met, which
    # are frequently consumed downstream.
    kwargs['start_met'] = string_to_met_utc(kwargs.get('startdate'), lazy=True)
    kwargs['stop_met'] = kwargs.get('start_met') + kwargs.get('duration')
    # Build the observation timeline.
    kwargs['timeline'], kwargs['gti_list'], calib_pattern, = _build_timeline(roi_model, **kwargs)
    kwargs['calib_pattern'] = calib_pattern
    return outfile, random_seed, roi_model, chrono, calib_pattern


def xpobssim(**kwargs):
    """Run the ixpeobssim fast simulator.
    """
    outfile, random_seed, roi_model, chrono, calib_pattern = _prepare_simulation(kwargs)
    # Setup the charging.
    _update_charging_settings(kwargs)
    # Loop over the detector units.
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
        irf_set = load_irf_set(kwargs.get('irfname'), du_id, gray_filter=kwargs.get('grayfilter'))
        logger.info('Done %s.', chrono)
        # Run the actual simulation.
        logger.info('Generating the photon list...')
        event_list = roi_model.rvs_event_list(irf_set, **kwargs)
        # Retrieve the onboard calibration intervals for the DU being simulated.
        calib_runs = calib_pattern[du_id]
        if len(calib_runs) > 0:
            logger.info('Activating calibration sources for DU %d...', du_id)
            calib_event_list = xEventList()
            # Loop over the onboard calibration intervals.
            for run in calib_runs:
                logger.info('Starting calibration @ MET %.3f--%.3f (%.3f s)...',
                            run.start_met, run.stop_met, run.duration)
                _kwargs = dict(start_met=run.start_met, duration=run.duration,
                               deadtime=kwargs.get('deadtime'))
                calib_event_list += run.calibration_source.rvs_event_list(irf_set, **_kwargs)
            # If we have calibration events, we need to concatenate the event
            # list, sort them in time and recalculate the trigger_id column.
            if calib_event_list.num_events() > 0:
                event_list += calib_event_list
        event_list.write_fits('xpobssim.py', roi_model, irf_set, **kwargs)
        logger.info('Done for detector unit # %d %s.', du_id, chrono)
    logger.info('All done %s!', chrono)
    return outlist


def main():
    xpobssim(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()

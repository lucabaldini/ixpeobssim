#!/usr/bin/env python
#
# Copyright (C) 2015--2025, the ixpeobssim team.
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
"""Separate time windows when the detector is illuminated/not illuminated by
the sun.

This application uses housekeeping information to split observation files 
into an INSUN and an INECLIPSE section, based on sun visibility. The livetime 
of the two parts is updated based on LV1 data.

Required input folder structure
-------------------------------
The application expects a standard IXPE observation folder layout::

    <folder>/
        event_l2/   Level-2 calibrated event files
                    (matched by ``ixpe*_det<DU>_evt2_v??.fits*``)
        event_l1/   Level-1 event files used to recompute LIVETIME
                    (matched by ``ixpe*_det<DU>_evt1_v??.fits*``)
        hk/         Housekeeping files containing spacecraft attitude data
                    (matched by ``ixpe*_all_adc_0110_v??.fits*``)

xpsun supports both .fits and .fits.gz files, but the latter will prevent
memory mapping potentially leading to several tens of gigabites of RAM usage.
Consider uncompressing files for large data sets.

level-2 files can be passed explicitly via ``--l2files`` instead of using
automatic parsing; the level-1 and HK files are always located by scanning
the observation folder. Beware: level 1 files are all those contained in the
event_l1 folder, so if some processing products are present they will all
be included leading to potential issues. A clean folder structure is recommended
when running xpsun. In case of level 2 subselections, use --l2files to 
specify the files to process.

GTI creation
------------
The core of the sun/eclipse separation is the ``create_ineclipse_gtis()``
function (in ``ixpeobssim.instrument.sc``), which reads the HK file columns:

* **ADSEC2SUN** -- angular distance of the sun relative to the spacecraft
* **ADSEC2ECL** -- angular distance encoding the Earth-occultation geometry

The INECLIPSE condition is defined as:

    ``(ADSEC2SUN < 0) AND (ADSEC2ECL >= 0)``

while the INSUN condition is its complement:

    ``(ADSEC2ECL < 0) OR (ADSEC2SUN >= 0)``

The boolean mask over the HK TIME column is converted into arrays of GTI
start/stop pairs.

Event filtering
---------------
For each detector unit and each level-2 file, the ``intersect_gti()`` function
(in ``ixpeobssim.evt.event``) performs the following steps:

1. Reads the existing GTI table from the level-2 event file.
2. Intersects it with the INSUN (or INECLIPSE) GTI intervals, producing a
   merged set of good-time intervals.
3. Replaces the GTI extension in the FITS file with the new intervals.
4. Filters the event list, keeping only rows whose TIME falls within the
   new GTIs.
5. Recomputes LIVETIME by summing the LIVETIME column from matching level-1
   files over the new GTI intervals, and updates the ONTIME keyword from the
   total good time.
6. Writes the result to a new file with a ``_insun`` or ``_inecl`` tag
   appended to the original filename.

Output files
------------
For an input file ``ixpe01234567_det1_evt2_v02.fits``, the application
produces:

* ``ixpe01234567_det1_evt2_v02_insun.fits``  -- events in sun-illuminated GTIs
* ``ixpe01234567_det1_evt2_v02_inecl.fits``  -- events in eclipse GTIs

Both files retain the same FITS structure as the original level-2 file,
with updated GTI, LIVETIME and ONTIME values.
"""

import glob
import os

from ixpeobssim.binning.misc import xBinnedLightCurve
from ixpeobssim.core import pipeline
from ixpeobssim.evt.event import xEventFile, intersect_gti
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.instrument.sc import create_ineclipse_gtis
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.utils.os_ import filter_input_file_list


PARSER = xArgumentParser(description=__description__)
PARSER.add_argument('folder', type=str,
                    help='path to the input observation folder')
PARSER.add_argument('--l2files', type=str, default=None, nargs='+',
                    help='level 2 file list')
PARSER.add_argument('--show', action='store_true', default=False,
                    help='show INSUN/INECLIPSE bands over the observation '\
                         'light curve')
PARSER.add_overwrite()
       

def build_l2_file_dict_from_folder(folder):
    """ Note: this is sorting the files alphabetically, in case this is
    relevant in any way
    """
    l2_folder = os.path.join(folder, 'event_l2')
    logger.debug('Searching for LEVEL 2 files in %s' % l2_folder)
    file_dict = {}
    for du_id in DU_IDS:
        match_list = glob.glob(
            os.path.join(l2_folder, 'ixpe*_det%s_evt2_v??.fits*' % du_id)
        )
        file_list = sorted(filter_input_file_list(match_list))
        logger.debug('%d LV2 files found for DU %d: %s' \
                     % (len(file_list), du_id, file_list))
        file_dict[du_id] = file_list
    return file_dict


def build_l2_file_dict_from_file_list(file_list):
    """ Note: this preserve the order of the files given by the user, in case
    this is relevant in any way
    """
    file_dict = {}
    for file_path in file_list:
        if 'det' not in file_path:
            raise ValueError('Cannot find a valid DU number for file %s' % file_path)
        file_name = os.path.basename(file_path)
        pieces = file_name.split('det')
        print(pieces[1])
        du_id = int(pieces[1][0])
        if (du_id not in DU_IDS):
            raise ValueError('Cannot find a valid DU number for file %s' % file_path)
        if du_id not in file_dict:
            file_dict[du_id] = [file_path]
        else:
            file_dict[du_id].append(file_path)
    return file_dict


def build_l1_file_dict_from_folder(folder, du_ids=DU_IDS):
    """ Note: this is sorting the files alphabetically, in case this is
    relevant in any way
    """
    l1_folder = os.path.join(folder, 'event_l1')
    logger.debug('Searching for LEVEL 1 files in %s' % l1_folder)
    file_dict = {}
    for du_id in du_ids:
        match_list = glob.glob(
            os.path.join(l1_folder, 'ixpe*_det%s_evt1_v??.fits*' % du_id)
        )
        file_list = sorted(filter_input_file_list(match_list))
        logger.debug('%d LV1 files found for DU %d: %s' \
                     % (len(file_list), du_id, file_list))
        file_dict[du_id] = file_list
    return file_dict


def find_hk_file_list_in_folder(folder):
    """
    """
    hk_folder = os.path.join(folder, 'hk')
    logger.debug('Searching for HK files in %s' % hk_folder)
    match_list = glob.glob(
        os.path.join(hk_folder, 'ixpe*_all_adc_0110_v??.fits*')
    )
    file_list = sorted(filter_input_file_list(match_list))
    if len(file_list) == 0:
        raise ValueError('Cannot find the appropriate hk file(s) in '\
                         'folder %s' % hk_folder)
    logger.debug('%d hk files found: %s'  % (len(file_list), file_list))
    return file_list


def make_plot(obs_path, inecl_path, insun_path, figname=''):
    """
    """
    inecl_file = xEventFile(inecl_path)
    insun_file = xEventFile(insun_path)
    plt.figure(f'{figname}')
    plot_gtis(inecl_file.gti_data(), label='INECLIPSE', hatch='/',
              edgecolor=None)
    plot_gtis(insun_file.gti_data(), color='r', label='INSUN', hatch='x',
              edgecolor=None)
    light_curve_path = pipeline.xpbin(obs_path, tbins=100000, 
                                      algorithm='LC', overwrite=True)
    light_curve = xBinnedLightCurve.from_file_list(light_curve_path)
    rate, rate_error = light_curve.rate(), light_curve.rate_error()
    plt.ylabel('Rate [Hz]')
    plt.xlabel('MET[s]')
    plt.errorbar(light_curve.TIME, rate, rate_error, fmt='o')


def plot_gtis(gti_data, color='g', alpha=0.3, label=None, **plot_opts):
    """
    """
    starts = gti_data['START']
    stops = gti_data['STOP']
    for start, stop in zip(starts, stops):
        span = plt.axvspan(start, stop, color=color, alpha=alpha, **plot_opts)
    if label is not None:
        span.set_label(label)
        plt.legend()

def xpsun(**kwargs):
    """ xpsun main app splitting the level 2 file into insun and ineclipse parts, 
    optionally plotting the light curve with the GTIs overlaid.
    """
    folder = kwargs['folder']
    show = kwargs.get('show', False)
    l2files = kwargs.get('l2files', None)
    overwrite=kwargs.get('overwrite', False)
    if l2files is not None:
        l2_file_dict = build_l2_file_dict_from_file_list(l2files)
    else:
        l2_file_dict = build_l2_file_dict_from_folder(folder)
    
    l1_file_dict = build_l1_file_dict_from_folder(folder,
                                                  du_ids=l2_file_dict.keys())
    hk_file_list = find_hk_file_list_in_folder(folder)

    inecl_paths = []
    insun_paths = []

    for du_id, l2_file_paths in l2_file_dict.items():
        l1_file_paths = l1_file_dict[du_id]
        for l2_file_path in l2_file_paths:
            logger.info('Processing file %s ...' % l2_file_path)
            base = os.path.splitext(l2_file_path)[0]
            # These are necessary for overwrite and are redefined in intersect_gti, but we need 
            # them here to check for existing files before running the processing
            inecl_path = base + '_inecl.fits'
            insun_path = base + '_insun.fits'
            if os.path.exists(inecl_path) and not overwrite:
                logger.warning('Output file %s already exists, skipping.', inecl_path)
            else:
                inecl_starts, inecl_stops = create_ineclipse_gtis(*hk_file_list)
                inecl_path = intersect_gti(l2_file_path, l1_file_paths, 
                                       inecl_starts, inecl_stops,
                                       tag='inecl')
            inecl_paths.append(inecl_path)
            if os.path.exists(insun_path) and not overwrite:
                logger.warning('Output file %s already exists, skipping.', insun_path)
            else:
                insun_starts, insun_stops = create_ineclipse_gtis(*hk_file_list,
                                                                  complement=True)
                insun_path = intersect_gti(l2_file_path, l1_file_paths, 
                                           insun_starts, insun_stops,
                                           tag='insun')
            insun_paths.append(insun_path)
            if show:    
                make_plot(l2_file_path, inecl_path, insun_path,
                          figname=f'DU{du_id}, {l2_file_path}')
    if show:
        plt.show()
    return inecl_paths, insun_paths

def main():
    """
    """
    args = PARSER.parse_args()
    xpsun(**args.__dict__)


if __name__ == '__main__':
    main()    

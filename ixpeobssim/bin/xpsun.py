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


PARSER = xArgumentParser()
PARSER.add_argument('folder', type=str,
                    help='path to the input observation folder')
PARSER.add_argument('--l2files', type=str, default=None, nargs='+',
                    help='level 2 file list')
PARSER.add_argument('--show', action='store_true', default=False,
                    help='show INSUN/INECLIPSE bands over the observation '\
                         'light curve')
       

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


def main():
    """
    """
    args = PARSER.parse_args()
    show = args.show
    if args.l2files is not None:
        l2_file_dict = build_l2_file_dict_from_file_list(args.l2files)
    else:
        l2_file_dict = build_l2_file_dict_from_folder(args.folder)
    
    l1_file_dict = build_l1_file_dict_from_folder(args.folder,
                                                  du_ids=l2_file_dict.keys())
    hk_file_list = find_hk_file_list_in_folder(args.folder)
    for du_id, l2_file_paths in l2_file_dict.items():
        l1_file_paths = l1_file_dict[du_id]
        for l2_file_path in l2_file_paths:
            logger.info('Processing file %s ...' % l2_file_path)
            inecl_starts, inecl_stops = create_ineclipse_gtis(*hk_file_list)
            inecl_path = intersect_gti(l2_file_path, l1_file_paths, 
                                       inecl_starts, inecl_stops,
                                       tag='inecl')
            
            
            insun_starts, insun_stops = create_ineclipse_gtis(*hk_file_list,
                                                              complement=True)
            insun_path = intersect_gti(l2_file_path, l1_file_paths, 
                                       insun_starts, insun_stops,
                                       tag='insun')
            if show:    
                make_plot(l2_file_path, inecl_path, insun_path,
                          figname=f'DU{du_id}, {l2_file_path}')
    if show:
        plt.show()


if __name__ == '__main__':
    main()    

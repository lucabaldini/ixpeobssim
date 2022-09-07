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


__description__ = \
"""Extract charging maps from an observation file.

This application is meant to open an observation file, read the CHRG_MAP
extension and copy it in a different FITS file, adding the proper header to
the PRIMARY extension so that the format is the same as the actual charging map
files.
"""

import os

from astropy.io import fits

from ixpeobssim.utils.time_ import current_datetime_string_utc, string_to_met_utc, met_to_string_utc
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.instrument.charging import xChargingPrimaryHDU
from ixpeobssim.utils.argparse_ import xArgumentParser


def extract_charging_map_file(input_file, outfolder=None, overwrite=True):
    """Read the input event file and extract the charging map.
    """
    hdu_list = fits.open(input_file)
    primary_header = hdu_list['PRIMARY'].header
    charging_hdu = hdu_list['CHRG_MAP']
    detnam = primary_header['DETNAM']
    det_id = primary_header['DET_ID']
    keywords = {
        'CREATOR' : 'xpchrgmap',
        'DATE'    : current_datetime_string_utc(fmt='%m/%d/%Y %H:%M:%S'),
        'DETNAM'  : detnam,
        'DET_ID'  : det_id
    }
    primary_hdu = xChargingPrimaryHDU()
    for key, value in keywords.items():
        primary_hdu.set_keyword(key, value)
    charging_hdu.header['DETNAM'] = detnam
    charging_hdu.header['DET_ID'] = det_id
    new_hdul = fits.HDUList([primary_hdu, charging_hdu])
    base_name = os.path.basename(input_file)
    if outfolder is None:
        outfolder = os.path.dirname(input_file)
    map_date = string_to_met_utc(charging_hdu.header['CVSD0001'], fmt='%m/%d/%Y')
    map_date_str = met_to_string_utc(map_date, fmt='%Y%m%d')
    output_file_name = 'ixpe_%s_chrgmap_%s' % (map_date_str, base_name)
    output_file_path = os.path.join(outfolder, output_file_name)
    logger.info('Writing charging map to %s...', output_file_path)
    new_hdul.writeto(output_file_path, overwrite=overwrite)
    logger.info('Done.')
    return output_file_path


def xpchrgmap(**kwargs):
    """Worker function.
    """
    file_list = kwargs.pop('filelist')
    return [extract_charging_map_file(file_path, **kwargs) for file_path in file_list]



"""Command-line switches."""
PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_outfolder()
PARSER.add_overwrite()



def main():
    xpchrgmap(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

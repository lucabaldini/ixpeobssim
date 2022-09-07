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

"""Conversion facilities for the legacy response files.

This was added after https://bitbucket.org/ixpesw/ixpeobssim/issues/564/
and serves the only purpose of providing facilities to update old response
files to the new formats without changing any essential content.

A few notes about hydiosynchrasies that were noticed through the release process:

* ixpemc_nocut_v1 only came for DU 1---we created copies for DU 2 and 3;
* the first three iterations of the response functions came with no modulation
  response functions, and we created the latter from the effective area and the
  modulation factor;
* the first three iterations of the response files shipped effective area curves
  combined for the three DU, and for the mirrors---those are essentially
  useless duplicated that served no purpose and were deleted.

For completeness, here is the full list of reponse files that were deleted in
the process:

.. code-block::

   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu123mmav002.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu123mmav003.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu123stdcutv001.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu123stdcutv002.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu123stdcutv003.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu1mmav001.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu1mmav002.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu1mmav003.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu2mmav002.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu2mmav003.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu3mmav002.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/arf/ixpemcdu3mmav003.arf'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/modfact/ixpemcdu123modfstdcutv001.fits'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/modfact/ixpemcdu123modfstdcutv002.fits'
   rm 'ixpeobssim/caldb/ixpe/gpd/cpf/modfact/ixpemcdu123modfstdcutv003.fits'
   rm 'ixpeobssim/caldb/ixpe/xrt/bcf/vign/ixpemcdu123vignstdcutv001.fits'
   rm 'ixpeobssim/caldb/ixpe/xrt/bcf/vign/ixpemcdu123vignstdcutv002.fits'
   rm 'ixpeobssim/caldb/ixpe/xrt/bcf/vign/ixpemcdu123vignstdcutv003.fits'
"""

from __future__ import print_function, division

import os

from ixpeobssim.core.fitsio import read_hdu_list_in_memory
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.instrument.du import du_logical_name, du_physical_name
from ixpeobssim.irf.caldb import IRF_TYPES, irf_file_path
from ixpeobssim.utils.logging_ import logger


# pylint: disable=no-member


_LEGACY_IRF_NAME_DICT = {
    'ixpemc_nocut_v1'   : 'ixpe:legacy:v1',
    'ixpemc_stdcut_v1'  : 'ixpe:legacy_stdcut:v1',
    'ixpemc_nocut_v2'   : 'ixpe:legacy:v2',
    'ixpemc_stdcut_v2'  : 'ixpe:legacy_stdcut:v2',
    'ixpemc_nocut_v3'   : 'ixpe:legacy:v3',
    'ixpemc_stdcut_v3'  : 'ixpe:legacy_stdcut:v3',
    'ixpemc_stdcut_v4'  : 'ixpe:legacy_stdcut:v4',
    'ixpemc_stdcut_v5'  : 'ixpe:legacy_stdcut:v5',
    'ixpemc_stdcut_v6'  : 'ixpe:legacy_stdcut:v6',
    'ixpe_obssim_v7'    : 'ixpe:legacy:v7',
    'ixpe_mc_v8'        : 'ixpe:legacy:v8',
    'ixpe_mcalpha075_v8': 'ixpe:legacy_alpha075:v8',
    'ixpe_mc_v9'        : 'ixpe:legacy:v9',
    'ixpe_mcalpha075_v9': 'ixpe:legacy_alpha075:v9'
}



def check_response_files():
    """Simple loop over the response files to identify possible orphans.
    """
    for irf_name in _LEGACY_IRF_NAME_DICT.values():
        for irf_type in IRF_TYPES:
            for du_id in DU_IDS:
                file_path = irf_file_path(irf_name, du_id, irf_type, check_file=False)
                if not os.path.exists(file_path):
                    logger.warning('Cannot find %s', file_path)


def _create_response_copy(irf_name, irf_type, src_du_id, dest_du_id):
    """Create a copy of a given response file for a given DU into a new file for
    a different DU.
    """
    assert src_du_id != dest_du_id
    src = irf_file_path(irf_name, src_du_id, irf_type, check_file=True)
    dest = irf_file_path(irf_name, dest_du_id, irf_type, check_file=False)
    if os.path.exists(dest):
        logger.warning('Output file %s exists, skipping...', dest)
        return
    hdu_list = read_hdu_list_in_memory(src, extension=None)
    logger.info('Setting the DU header keywords...')
    for hdu in hdu_list:
        hdu.header['DETNAM'] = du_logical_name(dest_du_id)
        hdu.header['DET_ID'] = du_physical_name(dest_du_id)
    logger.info('Writing modified HDU list to %s...', dest)
    hdu_list.writeto(dest)


def _create_rmf_files(irf_name):
    """Create a modulation response file from the effective area and modulation factor.
    """
    for du_id in DU_IDS:
        dest = irf_file_path(irf_name, du_id, 'mrf', check_file=False)
        if os.path.exists(dest):
            logger.warning('Output file %s exists, skipping...', dest)
            continue
        hdu_list = read_hdu_list_in_memory(irf_file_path(irf_name, du_id, 'arf'), 'arf')
        modf_hdu_list = read_hdu_list_in_memory(irf_file_path(irf_name, du_id, 'modf'))
        logger.info('Updating the response file...')
        for hdu in hdu_list:
            hdu.header['IRFTYPE'] = 'mrf'
        hdu_list[1].data['SPECRESP'] = hdu_list[1].data['SPECRESP'] * \
            modf_hdu_list[1].data['MODFRESP']
        hdu_list.writeto(dest)


def fill_missing_response_files():
    """Create the missing response files.

    This amounts to the response files for DU 2 and 3 for the ixpe:legacy:v1
    iteration and the modulation response files for the first three iterations
    of the response files.
    """
    for irf_name in ('ixpe:legacy:v1',):
        for irf_type in ('arf', 'rmf', 'modf', 'vign', 'psf'):
            _create_response_copy(irf_name, irf_type, 1, 2)
            _create_response_copy(irf_name, irf_type, 1, 3)
    for irf_name in ('ixpe:legacy:v1', 'ixpe:legacy_stdcut:v1', 'ixpe:legacy:v2', 'ixpe:legacy:v3'):
        _create_rmf_files(irf_name)


def format_response_files():
    """Loop over all the old response files (after they have been renamed according
    to the new naming conventions) and fix all the inconsistencies in the
    internal structures---strictly preserving the data.
    """
    for irf_name in _LEGACY_IRF_NAME_DICT.values():
        for irf_type in IRF_TYPES:
            for du_id in DU_IDS:
                file_path = irf_file_path(irf_name, du_id, irf_type, check_file=False)
                if not os.path.exists(file_path):
                    logger.warning('Cannot find %s', file_path)
                    continue
                hdu_list = read_hdu_list_in_memory(file_path, extension=None)
                hdu_list_modified = False

                # Check the file name.
                #for hdu in hdu_list:
                #    print(hdu.header.get('FILENAME'))

                # Check the DETNAM and DET_ID header keywords---note these have to
                # be in all the extensions.
                for hdu in hdu_list:
                    try:
                        assert hdu.header['DETNAM'] == du_logical_name(du_id)
                    except KeyError:
                        logger.warning('DETNAM not found in %s...', file_path)
                        hdu.header['DETNAM'] = du_logical_name(du_id)
                        hdu_list_modified = True
                    try:
                        assert hdu.header['DET_ID'] == du_physical_name(du_id)
                    except KeyError:
                        logger.warning('DET_ID not found in %s...', file_path)
                        hdu.header['DET_ID'] = du_physical_name(du_id)
                        hdu_list_modified = True

                # Check for old-style extention name in modf files---this only
                # applies to the extension at index 1.
                if hdu_list[1].header['EXTNAME'] == 'MODFRESP':
                    logger.warning('Old MODFRESP extension name for %s...', file_path)
                    hdu_list[1].header['EXTNAME'] = 'SPECRESP'
                    hdu_list_modified = True

                # Write the output file.
                if hdu_list_modified:
                    logger.info('Writing modified reponse file to %s...', file_path)
                    hdu_list.writeto(file_path, overwrite=True)



if __name__ == '__main__':
    #check_response_files()
    format_response_files()

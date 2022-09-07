# Copyright (C) 2016--2022, the ixpeobssim team.
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

"""Basic IRF naming conventions and CALDB-related facilities.
"""

from __future__ import print_function, division

import os

from ixpeobssim import IXPEOBSSIM_CALDB
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.os_ import check_input_file


# pylint: disable=invalid-name


# Basic CALDB folder structure.
__CALDB_FOLDER_DICT = {
    'arf' : ('gpd', 'cpf', 'arf'),
    'rmf' : ('gpd', 'cpf', 'rmf'),
    'mrf' : ('gpd', 'cpf', 'mrf'),
    'modf': ('gpd', 'cpf', 'modfact'),
    'psf' : ('xrt', 'bcf', 'psf'),
    'vign': ('xrt', 'bcf', 'vign')
    }


# And list of all the available IRF types.
IRF_TYPES = __CALDB_FOLDER_DICT.keys()



def irf_file_name(base, du_id, irf_type, intent, version):
    """Return the file name for a specific response function.

    Arguments
    ---------
    base : str
        Base name of the set of response functions.

    du_id : int
        Identifier of the specific DU.

    irf_type : str
        Identifier for the calibration data.

    intent : str
        The specific intent for the response file.

    version : int
        Version number.

    The basic naming convention scheme used for the CALDB by HEASARC is
    ixpe_[instrument]_[datatype]_[date]_[ver].[ext], where:

    * [instrument] indicates the detector (d1/d2/d3);
    * [datatype] provides an identifier for the calibration data;
    * [date] indicates the first date of validity of the file;
    * [ver] is the CALDB version number for that file.

    We approximately follow the HEASARC conventions, except for the fact that
    we don't really have a concept of validity date.
    """
    if irf_type in ('arf', 'mrf', 'rmf'):
        return '%s_d%d_%s_v%03d.%s' % (base, du_id, intent, version, irf_type)
    if irf_type == 'modf':
        irf_type = 'mfact'
    return '%s_d%d_%s_%s_v%03d.fits' % (base, du_id, intent, irf_type, version)


def parse_irf_name(irf_name, delimiter=':'):
    """Parse a generic IRF name and return the basic field values it encapsulates.

    The ``irf_name`` is an internal designation that ixpeobssim uses to
    identitify a full (self-consistent) set of response functions.
    """
    try:
        base, intent, version = irf_name.split(delimiter)
    except ValueError as e:
        logger.error('Error in parse_irf_name(): %s', e)
        abort('Invalid IRF name %s' % irf_name)
    version = int(version.strip('v'))
    return base, intent, version


def irf_folder_path(irf_type, caldb_path=IXPEOBSSIM_CALDB):
    """Return the CALDB folder for a particular IRF type.
    """
    if not irf_type in __CALDB_FOLDER_DICT:
        abort('Invalid IRF type (%s)' % irf_type)
    return os.path.join(caldb_path, 'ixpe', *__CALDB_FOLDER_DICT[irf_type])


def irf_file_path(irf_name, du_id, irf_type, caldb_path=None, check_file=True):
    """Return the full file path to a particular IRF file.
    """
    if caldb_path is None:
        caldb_path = IXPEOBSSIM_CALDB
    folder_path = irf_folder_path(irf_type, caldb_path)
    base, intent, version = parse_irf_name(irf_name)
    file_name = irf_file_name(base, du_id, irf_type, intent, version)
    file_path = os.path.join(folder_path, file_name)
    if check_file:
        check_input_file(file_path)
    return file_path

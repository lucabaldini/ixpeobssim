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

"""IRF top-level utilities.

This module contains all the functions to load response function, as well as
a convenience class, xIRFSet, encapsulating a coherent set of files and the
corresponding method to load the complete set.

In addition, the peek_irf_type() convenience method allows to check the type
of the response file, based on the IRFTYPE keyword of the primary header.

Note that, by default, all the response objects are cached in memory, so that the
corresponding files are not read from disk multiple times within the span of a
single program---this speeds up things when the same operation (e.g., the binning
of an event file) is repeated many times, such as in benchmarks.
"""

from __future__ import print_function, division


import os

from astropy.io import fits
import numpy

from ixpeobssim import IXPEOBSSIM_CALDB
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.irf.arf import xEffectiveArea, xTowEffectiveArea
from ixpeobssim.irf.caldb import irf_file_path, parse_irf_name
from ixpeobssim.irf.legacy import _LEGACY_IRF_NAME_DICT
from ixpeobssim.irf.modf import xModulationFactor
from ixpeobssim.irf.mrf import xModulationResponse
from ixpeobssim.irf.psf import xPointSpreadFunction
from ixpeobssim.irf.rmf import xEnergyDispersion
from ixpeobssim.irf.vign import xVignetting
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.os_ import check_input_file


# pylint: disable=invalid-name, no-member, too-many-arguments, too-many-instance-attributes


# Name of the IRF set to be used by default throughout the package.
DEFAULT_IRF_NAME = 'ixpe:obssim:v12'

# Private dictionary to cache the objects that have already been loaded.
__CACHE = {}


# pylint: disable=inconsistent-return-statements
def peek_irf_type(file_path):
    """Open a given response file and retrieve the IRF type it contains.
    """
    check_input_file(file_path)
    with fits.open(file_path) as hdu_list:
        try:
            return hdu_list['PRIMARY'].header['IRFTYPE']
        except KeyError:
            abort('IRFTYPE keyword not found in header of %s.' % (file_path))


def _load_irf_base(cls, irf_type, irf_name=DEFAULT_IRF_NAME, du_id=1,
    caldb_path=None, cache=True, simple_weighting=False, gray_filter=False):
    """Base helper function to load a response function.

    Arguments
    ---------
    cls : class
        The class of the response function to be loaded.

    irf_type : str
        The type of response file.

    irf_name : str
        The name of the response file.

    du_id : int
        The logical designator of the detector unit.

    caldb_path : str
        The path to the calibration database.

    cache : bool
        If True, the response object is cached in memory for future use once loaded.

    simple_weighting : bool
        If True, load the response file with the SIMPLE weighting scheme.

    gray_filter : bool
        If True, load the response files apppropriate for the gray filtes (where)
        this makes sense.
    """
    # Small hook to support old-style IRF names.
    if irf_name in _LEGACY_IRF_NAME_DICT:
        _irf_name = irf_name
        irf_name = _LEGACY_IRF_NAME_DICT[_irf_name]
        logger.info('Old-style IRF name %s -> %s', _irf_name, irf_name)
    file_path = irf_file_path(irf_name, du_id, irf_type, caldb_path, True,
        simple_weighting, gray_filter)
    if file_path in __CACHE:
        logger.info('Using cached %s object at %s...', cls.__name__, file_path)
        return __CACHE[file_path]
    irf = cls(file_path)
    if cache:
        __CACHE[file_path] = irf
    return irf


def load_arf(irf_name=DEFAULT_IRF_NAME, du_id=1, caldb_path=None, cache=True,
    simple_weighting=False, gray_filter=False):
    """Facility to load the effective area for a given IRF set.
    """
    return _load_irf_base(xEffectiveArea, 'arf', irf_name, du_id, caldb_path,
        cache, simple_weighting, gray_filter)


def load_vign(irf_name=DEFAULT_IRF_NAME, du_id=1, caldb_path=None, cache=True):
    """Facility to load the vignetting for a given IRF set.
    """
    return _load_irf_base(xVignetting, 'vign', irf_name, du_id, caldb_path, cache)


def load_psf(irf_name=DEFAULT_IRF_NAME, du_id=1, caldb_path=None, cache=True):
    """Facility to load the point-spread function for a given IRF set.
    """
    return _load_irf_base(xPointSpreadFunction, 'psf', irf_name, du_id, caldb_path, cache)


def load_modf(irf_name=DEFAULT_IRF_NAME, du_id=1, caldb_path=None, cache=True):
    """Facility to load the modulation factor for a given IRF set.
    """
    return _load_irf_base(xModulationFactor, 'modf', irf_name, du_id, caldb_path, cache)


def load_mrf(irf_name=DEFAULT_IRF_NAME, du_id=1, caldb_path=None, cache=True,
    simple_weighting=False, gray_filter=False):
    """Facility to load the modulation response for a given IRF set.
    """
    return _load_irf_base(xModulationResponse, 'mrf', irf_name, du_id, caldb_path,
        cache, simple_weighting, gray_filter)


def load_rmf(irf_name=DEFAULT_IRF_NAME, du_id=1, caldb_path=None, cache=True):
    """Facility to load the energy dispersion for a given IRF set.
    """
    return _load_irf_base(xEnergyDispersion, 'rmf', irf_name, du_id, caldb_path, cache)


def load_tow_aeff(irf_name=DEFAULT_IRF_NAME, du_id=1, caldb_path=None, cache=True):
    """Facility to load the effective area for a given IRF set.
    """
    return _load_irf_base(xTowEffectiveArea, 'tow', irf_name, du_id, caldb_path, cache)



class xIRFSet:

    """Small container class representing a full set of IRFs.
    """

    # pylint: disable=too-few-public-methods

    def __init__(self, irf_name, du_id, caldb_path=None, cache=True,
        simple_weighting=False, gray_filter=False):
        """Constructor.
        """
        self.irf_name = irf_name
        self.du_id = du_id
        args = irf_name, du_id, caldb_path, cache
        self.aeff = load_arf(*args, simple_weighting, gray_filter)
        self.vign = load_vign(*args)
        self.edisp = load_rmf(*args)
        self.psf = load_psf(*args)
        self.modf = load_modf(*args)
        self.mrf = load_mrf(*args, simple_weighting, gray_filter)



def load_irf_set(irf_name=DEFAULT_IRF_NAME, du_id=1, caldb_path=None, cache=True,
    simple_weighting=False, gray_filter=False):
    """Load a set of instrument response functions.
    """
    return xIRFSet(irf_name, du_id, caldb_path, cache, simple_weighting, gray_filter)

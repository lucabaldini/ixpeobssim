#!/usr/bin/env python
#
# Copyright (C) 2020, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public Licensese as published by
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

"""Common facilities for source models.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.system_ import import_module
from ixpeobssim.irf.ebounds import ENERGY_MIN, ENERGY_MAX



def load_tabular_data(file_path, emin=ENERGY_MIN, emax=ENERGY_MAX, energy_column=0, delimiter=None):
    """Generic routine to parse tabular data from ASCII files, see
    https://bitbucket.org/ixpesw/ixpeobssim/issues/279

    This is intended as a general purpose routine for all kind of models (e.g.,
    spectral and polarimetric), and it comes from the basic realization that
    we were re-implementing the same thing over and over again in the source
    configuration files. You can use this whenever you have a data file where
    one of the columns (typically the first) indicates the energy, and the
    others contain any kind of modeling information, e.g., the flux, or the
    polarization degree and angle.

    The method support a filtering between a minimum and maximum energy, which
    by default is performed between the energy bound of the response functions.

    Args
    ----
    file_path : str
        The path to the ASCII file containing the data.

    emin : float
        The minimum energy for filtering the data.

    emax : float
        The maximum energy for filtering the data.

    energy_column : int, default to 0
        The identifier of the column containing the energy values.
    """
    logger.info('Loading tabular data from %s...', file_path)
    data = numpy.loadtxt(file_path, delimiter=delimiter, unpack=True)
    energy = data[energy_column]
    logger.info('Done, %d columns and %d rows read out.', len(data), energy.size)
    mask = numpy.logical_and(energy >= emin, energy <= emax)
    return tuple([col[mask] for col in data])


def import_roi(file_path):
    """Import an ROI object from a python module.

    Arguments
    ---------
    file_path : string
        The path to the Python module.
    """
    return import_module(file_path).ROI_MODEL

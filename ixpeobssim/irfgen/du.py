#!/usr/bin/env python
#
# Copyright (C) 2019--2022, the ixpeobssim team.
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

from ixpeobssim import IXPEOBSSIM_IRFGEN
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.utils.logging_ import logger


IRFGEN_DU_DATA = os.path.join(IXPEOBSSIM_IRFGEN, 'data', 'du')

UVF_TRAN_FILE_NAME = 'uv_filter_nominal_transmission.txt'
UVF_TRAN_FILE_PATH = os.path.join(IRFGEN_DU_DATA, UVF_TRAN_FILE_NAME)


def uv_filter_transparency_spline():
    """Return a spline with the transparency of the UV filter as a function
    of the energy.
    """
    file_path = UVF_TRAN_FILE_PATH
    logger.info('Reading UV filter transparency data from %s...' % file_path)
    energy, trans = numpy.loadtxt(file_path, unpack=True)
    fmt = dict(xlabel='Energy [keV]', ylabel='UV filter transparency')
    return xInterpolatedUnivariateSplineLinear(energy, trans, **fmt)


def uv_filter_transparency(energy):
    """Return the transparency of the UV filter evaluated on a given grid of
    energy points.
    """
    trans = uv_filter_transparency_spline()(energy)
    assert (trans > 0).all()
    return trans

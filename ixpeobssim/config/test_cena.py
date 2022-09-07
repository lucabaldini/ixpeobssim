#!/usr/bin/env python
#
# Copyright (C) 2015--2016, the ixpeobssim team.
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

from ixpeobssim.srcmodel.roi import xChandraObservation, xChandraROIModel
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim import IXPEOBSSIM_TEST, IXPEOBSSIM_CONFIG_REG
from ixpeobssim.utils.astro import read_ds9


EVENT_FILE_PATH = os.path.join(IXPEOBSSIM_TEST, 'data', 'cena.fits')
REG_SOURCE_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG_REG, 'cena_jet+core.reg')
regions = read_ds9(REG_SOURCE_FILE_PATH)

ROI_MODEL = xChandraROIModel(EVENT_FILE_PATH, acis='I')

polarization_degree = constant(0.)
polarization_angle = constant(0.)

#Try to remove the core of the galaxy
core = xChandraObservation('Core', polarization_degree, polarization_angle,
                           regions[1], exclude=True)
ROI_MODEL.add_source(core)

cena = xChandraObservation('CenA', polarization_degree, polarization_angle)
ROI_MODEL.add_source(cena)

if __name__ == '__main__':
    print(ROI_MODEL)

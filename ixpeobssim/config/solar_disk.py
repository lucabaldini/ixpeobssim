#!/usr/bin/env python
#
# Copyright (C) 2015--2019, the ixpeobssim team.
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


import numpy, os

from ixpeobssim.srcmodel.roi import xUniformDisk, xROIModel
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.core.spline import xInterpolatedUnivariateLogSpline
from ixpeobssim.utils.logging_ import logger
from ixpeobssim import IXPEOBSSIM_CONFIG


RA = 10.
DEC = 10.
RAD = 16./60. #mean angular half size of the sun
#RAD = 0.5/60. #angular half size of a typical sunspot (active region)
POL_DEGREE = 0.
POL_ANGLE = numpy.radians(0.)

def parse_spectral_model(file_name, emin=0.5, emax=15.):
    """Parse the input file with the spectral point.
    """
    file_path = os.path.join(IXPEOBSSIM_CONFIG, 'ascii', file_name)
    logger.info('Parsing input file %s...' % file_path)
    energy, flux_flare, flux_quiet = numpy.loadtxt(file_path, unpack=True)
    _mask = (energy >= emin)*(energy <= emax)
    energy = energy[_mask]
    flux_quiet = flux_quiet[_mask]
    return energy, flux_quiet, flux_flare

energy, flux_quiet, flux_flare = parse_spectral_model('solar_spectrum.csv')

fmt = dict(xlabel='Energy [keV]', ylabel='Flux [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]')

quiet_solar_model = xInterpolatedUnivariateLogSpline(energy, flux_quiet,
                                                     **fmt)
flare_solar_model = xInterpolatedUnivariateLogSpline(energy, flux_flare,
                                                     **fmt)

def spec(E, t):
    return quiet_solar_model(E)
    #return flare_solar_model(E)

ROI_MODEL = xROIModel(RA, DEC)

polarization_degree = constant(POL_DEGREE)
polarization_angle = constant(POL_ANGLE)

solar_disk = xUniformDisk('Sun', RA, DEC, RAD, spec,
                          polarization_degree, polarization_angle)

ROI_MODEL.add_source(solar_disk)

def display():
    """Display the source model.
    """
    from ixpeobssim.utils.matplotlib_ import plt

    print(ROI_MODEL)
    plt.figure('Energy spectrum')
    quiet_solar_model.plot(logy=True, label='Quiet Sun')
    flare_solar_model.plot(logy=True, label='Flare')

    plt.legend(bbox_to_anchor=(0.35, 0.35))
    plt.show()

if __name__ == '__main__':
    display()

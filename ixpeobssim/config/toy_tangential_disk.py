#!/usr/bin/env python
#
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

"""
Simple uniform disk with a complex, tangential polarization pattern.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xUniformDisk, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import xTangentialPolarizationField
from ixpeobssim.utils.units_ import arcmin_to_degrees, degrees_to_arcmin
from ixpeobssim.utils.fmtaxis import fmtaxis


__model__ = file_path_to_model_name(__file__)

# Basic model parameters.
RA, DEC = 45., 45.
RADIUS = arcmin_to_degrees(5.)
PL_NORM = 1.
PL_INDEX = 2.
MAX_POL_DEG = 0.5

# Definition of the spectro-polarimetric properties
spec = power_law(PL_NORM, PL_INDEX)

def radial_profile(r, E, t=None):
    """Radial profile of the polarization degree.
    """
    emin = 0.
    emax = 8.
    return numpy.clip((MAX_POL_DEG * r / RADIUS) * (E - emin) / (emax - emin), 0., 1.)

pol_field = xTangentialPolarizationField(RA, DEC, radial_profile)
pol_deg = pol_field.polarization_degree_model()
pol_ang = pol_field.polarization_angle_model()

# Definition of the ROI.
src = xUniformDisk('Uniform disk', RA, DEC, RADIUS, spec, pol_deg, pol_ang)
ROI_MODEL = xROIModel(RA, DEC, src)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 100)
    # Energy spectrum
    plt.figure('%s spectrum' % __model__)
    plt.plot(energy, spec(energy))
    setup_gca(xmin=emin, xmax=emax, ymin=spec(emax), logx=True, logy=True, grids=True, **fmtaxis.spec)
    # Radial profile.
    r = numpy.linspace(0., RADIUS, 100)
    plt.figure('%s polarization radial profile' % __model__)
    for energy in (2., 4., 8.):
        plt.plot(degrees_to_arcmin(r), radial_profile(r, energy), label='%.1f keV' % energy)
    setup_gca(grids=True, legend=True, xlabel='Distance from center [arcmin]',
        ylabel='Polarization degree')



if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

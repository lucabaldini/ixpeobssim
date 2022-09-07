#!/usr/bin/env python
#
# Copyright (C) 2015--2022, the ixpeobssim team.
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


"""Full source model definition can be found in ixpeobssim.config.toy_casa.py.

The spectral model is taken from
E.A. Helder and J. Vink, "Characterizing the non-thermal emission of Cas A",
Astrophys.J. 686 (2008) 1094--1102, http://arxiv.org/abs/0806.3748. The spectrum
of Cas A is a complex superposition of thermal and non thermal emission, and for
our purposes, we call thermal anything that is making up for the lines and
non-thermal all the rest.

We have two images of Cas A, at low (1.5--3.0 keV) and high (4.0--6.0 keV) and,
due to the absence of lines between 4 and 6 keV we're attaching the latter to
the non-thermal spectrum and the former to the thermal component.
"""

from __future__ import print_function, division

import os

import numpy

from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII, IXPEOBSSIM_CONFIG_FITS
from ixpeobssim.config import bootstrap_display, file_path_to_model_name
from ixpeobssim.srcmodel.bkg import xTemplateInstrumentalBkg
from ixpeobssim.srcmodel.img import xFITSImage
from ixpeobssim.srcmodel.polarization import constant, xTangentialPolarizationField
from ixpeobssim.srcmodel.roi import xExtendedSource, xROIModel
from ixpeobssim.srcmodel.spectrum import load_spectral_spline
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.units_ import arcmin_to_degrees

# pylint: disable=invalid-name, unused-argument, consider-using-f-string

__model__ = file_path_to_model_name(__file__)


# Coordinates of the pointing.
RA_PNT = 350.8664167
DEC_PNT = 58.8117778
# Maximum polarization degree, and corresponding radius.
MAX_POL_DEG = 0.5
MAX_RADIUS = arcmin_to_degrees(2.8)

def _load_spec(file_name, emin=1., emax=15.):
    """Convenience function to load a spectral csv file.
    """
    file_path = os.path.join(IXPEOBSSIM_CONFIG_ASCII, file_name)
    return load_spectral_spline(file_path, emin=emin, emax=emax, delimiter=',')

# Read in the spectral models from the proper files.
total_spec_spline = _load_spec('casa_total_spectrum.csv')
non_therm_spec_spline = _load_spec('casa_nonthermal_spectrum.csv')
therm_spec_spline = total_spec_spline - non_therm_spec_spline

def therm_spec(E, t):
    """Definition of the photon spectrum for the thermal component.
    """
    return therm_spec_spline(E)

def non_therm_spec(E, t):
    """Definition of the photon spectrum for the non-thermal component.
    """
    return non_therm_spec_spline(E)

# The thermal component is unpolarized.
therm_pol_ang = constant(0.)
therm_pol_deg = constant(0.)

def non_therm_radial_profile(r, E=None, t=None):
    """Radial profile of the polarization degree.

    This is going from zero at the center of the field to the target maximum value
    at the specified radius. (Note the return value is clipped between 0 and 1
    to avoid unphysical values.)
    """
    return numpy.clip((MAX_POL_DEG * r / MAX_RADIUS), 0., 1.)

non_therm_pol_field = xTangentialPolarizationField(RA_PNT, DEC_PNT, non_therm_radial_profile)
non_therm_pol_deg = non_therm_pol_field.polarization_degree_model()
non_therm_pol_ang = non_therm_pol_field.polarization_angle_model()

# Finally, the Chandra images for the thermal and non thermal components.
le_img_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'casa_1p5_3p0_keV.fits')
he_img_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'casa_4p0_6p0_keV.fits')

# Create the model components.
therm_comp = xExtendedSource('Cas A thermal', le_img_file_path, therm_spec,
    therm_pol_deg, therm_pol_ang)
non_therm_comp = xExtendedSource('Cas A non-thermal', he_img_file_path, non_therm_spec,
    non_therm_pol_deg, non_therm_pol_ang)
instrumental_bkg = xTemplateInstrumentalBkg()

# Create the actual ROI object.
ROI_MODEL = xROIModel(RA_PNT, DEC_PNT)
ROI_MODEL.add_sources(therm_comp, non_therm_comp, instrumental_bkg)



def display():
    """Display the source model.
    """
    plt.figure('%s energy spectrum' % __model__)
    total_spec_spline.plot(label='Total')
    non_therm_spec_spline.plot(label='Non-thermal')
    therm_spec_spline.plot(label='Thermal')
    setup_gca(logy=True, legend=True, grids=True, xmin=1., xmax=15.)

    plt.figure('%s low-energy Chandra image' % __model__)
    therm_comp.image.plot(stretch='log')
    xFITSImage.add_label('Chandra 1.5-3.0 keV')

    plt.figure('%s hi-energy Chandra image' % __model__)
    non_therm_comp.image.plot(stretch='log')
    xFITSImage.add_label('Chandra 4.0-6.0 keV')



if __name__ == '__main__':
    bootstrap_display()

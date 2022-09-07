#!/usr/bin/env python
#
# Copyright (C) 2018--2019, the ixpeobssim team.
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

"""This is an example of a spectro-polarimetric model of a Compton-thick
AGN, and a well-known one, where we have detailed information about the source
geometry and energy spectrum.

There are three distinct spectral components to the model, each one with its
own polarimetric pattern:

* the cold reflection from the torus;
* the warm reflection from the cone;
* the emission lines.

.. image:: figures/models/ngc1068_energy_spectra.png

As far as the geometry for the polarimetric modeling is concerned, the 
ionization cone is assumed to be perpendicular to the plane of the torus, and
the inclination angle of the system is 75 degrees.

.. image:: figures/models/ngc1068_polarization_degree.png

We note that the position angles for the cold and warm relflection components
are orthogonal to each other, so that this source is a good illustration
of the harmonic addition in ixpeobssim.


Simulation output
~~~~~~~~~~~~~~~~~

Below is an output of an ixpeobssim simulation, where the measured polarization
degree, estimated in a few energy bins, is compared to the input model.

.. image:: figures/obssim/ngc1068_polarization_degree.png

References
~~~~~~~~~~

[1] Ghisellini, G.; Haardt, F.; Matt, G., "`The Contribution of the Obscuring Torus to the X-Ray Spectrum of Seyfert Galaxies - a Test for the Unification Model <http://adsabs.harvard.edu/abs/1994MNRAS.267..743G>`_", Monthly Notices of the Royal Astronomical Society, Vol. 267, NO. 3/APR1, P. 743, 1994

[2] Goosmann, R. W.; Matt, G., "`Modeling X-ray polarimetry while flying around the misaligned outflow of NGC 1068 <http://adsabs.harvard.edu/abs/2011sf2a.conf..583G>`_", SF2A-2011: Proceedings of the Annual meeting of the French Society of Astronomy and Astrophysics, 2011
"""

from __future__ import print_function, division

import os
import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.srcmodel.polarization import constant, constant_spline
from ixpeobssim.srcmodel.polarization import harmonic_component_addition


__model__ = file_path_to_model_name(__file__)
source_name = 'NGC 1068'
ra, dec = 40.6698792, -0.0132889


def parse_spectral_data(file_name, emin=1., emax=10.):
    """Parse the input data for the spectral components.
    
    The input format is a simple ascii file with two columns:
    * energy [keV]
    * flux [counts / s / cm2 /keV]
    """
    file_path = os.path.join(IXPEOBSSIM_CONFIG_ASCII, file_name)
    logger.info('Reading data from %s...' % file_path)
    energy, flux = numpy.loadtxt(file_path, unpack=True)
    mask = (energy >= emin) * (energy <= emax)
    energy = energy[mask]
    flux = flux[mask]
    fmt = dict(xlabel='Energy [keV]',
               ylabel='Flux [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]')
    return xInterpolatedUnivariateSplineLinear(energy, flux, **fmt)


def parse_polarization_data(file_name, emin=1., emax=10.):
    """Parse the input data for the polarization components.

    Again, the input format is a simple ascii file with two columns:
    * energy [keV]
    * polarization degree
    """
    file_path = os.path.join(IXPEOBSSIM_CONFIG_ASCII, file_name)
    logger.info('Reading data from %s...' % file_path)
    energy, pol_deg = numpy.loadtxt(file_path, unpack=True)
    mask = (energy >= emin) * (energy <= emax)
    energy = energy[mask]
    pol_deg = pol_deg[mask]
    fmt = dict(xlabel='Energy [keV]', ylabel='Polarization degree')
    return xInterpolatedUnivariateSpline(energy, pol_deg, k=3, **fmt)


# Parse all the necessary input data.
cold_spec_spline = parse_spectral_data('ngc1068_cold_spectrum.txt')
cold_pol_deg_spline = parse_polarization_data('ngc1068_cold_polarization.txt')
warm_spec_spline = parse_spectral_data('ngc1068_warm_spectrum.txt')
warm_pol_deg_spline = parse_polarization_data('ngc1068_warm_polarization.txt')
line_spec_spline = parse_spectral_data('ngc1068_lines_spectrum.txt')

# Text strings to label the components.
cold_label = 'Cold reflection'
warm_label = 'Warm reflection'
line_label = 'Lines'
total_label = '%s total' % source_name

# Define all the specral and polarimetric components: cold reflection...
def cold_spec(E, t):
    return cold_spec_spline(E)

def cold_pol_deg(E, t, ra, dec):
    return cold_pol_deg_spline(E)

cold_pol_ang = constant(numpy.radians(90.))

# ...warm reflection...
def warm_spec(E, t):
    return warm_spec_spline(E)

def warm_pol_deg(E, t, ra, dec):
    return warm_pol_deg_spline(E)

warm_pol_ang = constant(numpy.radians(0.))

# And, finally, lines...
def line_spec(E, t):
    return line_spec_spline(E)

line_pol_ang = constant(numpy.radians(0.))
line_pol_deg = constant(0.)

# Define a few more splines that, though not strictly necessary for the
# simulation of the event lists, might be useful.
energy = line_spec_spline.x
emin = energy[0]
emax = energy[-1]
deg_fmt = dict(ylabel='Polarization degree')
ang_fmt = dict(ylabel='Polarization angle [rad]')
spec_spline = cold_spec_spline + warm_spec_spline + line_spec_spline
cold_pol_ang_spline = constant_spline(emin, emax, cold_pol_ang(), **ang_fmt)
warm_pol_ang_spline = constant_spline(emin, emax, warm_pol_ang(), **ang_fmt)
line_pol_deg_spline = constant_spline(emin, emax, line_pol_deg(), **deg_fmt)
line_pol_ang_spline = constant_spline(emin, emax, line_pol_ang(), **ang_fmt)

# Harmonic addition of all the polarimetric components
cold_component = (cold_spec_spline, cold_pol_deg_spline, cold_pol_ang_spline)
warm_component = (warm_spec_spline, warm_pol_deg_spline, warm_pol_ang_spline)
line_component = (line_spec_spline, line_pol_deg_spline, line_pol_ang_spline)
components = (cold_component, warm_component, line_component)
_, pol_deg_spline, pol_ang_spline = harmonic_component_addition(*components)

# Move on creating the model components.
cold_src = xPointSource(cold_label, ra, dec, cold_spec, cold_pol_deg,
                        cold_pol_ang)
warm_src = xPointSource(warm_label, ra, dec, warm_spec, warm_pol_deg,
                        warm_pol_ang)
line_src = xPointSource(line_label, ra, dec, line_spec, line_pol_deg,
                        line_pol_ang)

# Create the complete region of interest.
ROI_MODEL = xROIModel(ra, dec, cold_src, warm_src, line_src)



def display(emin=1., emax=10.):
    """Display the source model.
    """
    # Energy spectrum.
    plt.figure('%s energy spectra' % __model__)
    cold_spec_spline.plot(label=cold_label)
    warm_spec_spline.plot(label=warm_label)
    line_spec_spline.plot(label=line_label)
    spec_spline.plot(label=total_label)
    setup_gca(logy=True, xmin=emin, xmax=emax, ymin=1e-6, legend=True)

    # Polarization degree as a function of the energy.
    plt.figure('%s polarization degree' % __model__)
    cold_pol_deg_spline.plot(label='%s (P.A. = %d$^\\circ$)' %\
                             (cold_label, numpy.degrees(cold_pol_ang())))
    warm_pol_deg_spline.plot(label='%s (P.A. = %d$^\\circ$)' %\
                             (warm_label, numpy.degrees(warm_pol_ang())))
    line_pol_deg_spline.plot(label=line_label)
    pol_deg_spline.plot(label=total_label)
    setup_gca(xmin=emin, xmax=emax, ymin=-0.01, legend=True)

    # Polarization angle as a function of the energy.
    plt.figure('%s polarization angle' % __model__)
    cold_pol_ang_spline.plot(label=cold_label)
    warm_pol_ang_spline.plot(label=warm_label)
    line_pol_ang_spline.plot(label=line_label)
    pol_ang_spline.plot(label=total_label)
    setup_gca(xmin=emin, xmax=emax, ymin=-0.01, legend=True)


if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

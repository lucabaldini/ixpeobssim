# Copyright (C) 2021, the ixpeobssim team.
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

"""

Spectro-polarimetric Characteristics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is plenty of useful information in
`Rea et al., 2003 <https://iopscience.iop.org/article/10.1086/374585>`_,
including the best fit parameters for the phase-averaged spectrum, i.e.,

* nH   : 1.36e22 cm^{-2}
* gamma: 2.4
* Flux : 1.87e-10 erg cm^{-1} s^{-1} in 0.5--10 keV

The basic spectral and polarimetric model comes from the dedicated model table,
courtesy of Roberto and Roberto. The parameters used for the simulations are


.. code-block::

    # Magnetar model parameters.
    MAG_MODEL = xMagnetarModelsT2020.BLACKBODY
    MAG_CHI = 85.
    MAG_XI = 55.
    MAG_DELTA_PHI = 0.35
    MAG_BETA = 0.39


Timing
~~~~~~

The main reference for the timing is
`Kuiper et al., 2006 <https://iopscience.iop.org/article/10.1086/504317>`_

The ephemeris for the source are given in table 3, and read:

* t0      : 52590.0 (MJD)
* nu0     : 0.09089812328 Hz
* nudot0  : -1.59836e-13 Hz s^-1
* nudddot : 0. Hz s^-2

"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name, bootstrap_display
from ixpeobssim.config.dc1_bkg import instrumental_bkg
from ixpeobssim.config.dc1_utils import align_ephemeris
from ixpeobssim.srcmodel.bkg import xExtragalacticBkg, xGalacticBkg
from ixpeobssim.srcmodel.magnetar import xMagnetarModelsT2020, xMagnetarTableModelT2020
from ixpeobssim.srcmodel.roi import xPeriodicPointSource, xROIModel
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color


#pylint: disable=invalid-name

__model__ = file_path_to_model_name(__file__)


# Source coordinates, in decimal degrees.
SRC_NAME = 'PSR J1708-4008'
SRC_RA, SRC_DEC = 257.20416667, -40.15277778
SRC_L, SRC_B = 346.48107483, 0.0277362

# Pointing coordinates
PNT_RA, PNT_DEC = SRC_RA, SRC_DEC

# Observation start date, and reference for the ephemeris.
START_DATE = '2022-09-01'
DURATION = 2000000.

# Interstellar absorption.
NH = 1.36e22

# Galactic diffuse.
GAL_BKG_RATE = 120.

# Phase-averaged integral flux.
INTEGRAL_FLUX = 1.87e-10
INT_FLUX_EMIN = 0.5
INT_FLUX_EMAX = 10

# Align the ephemeris to the start date.
EPHEMERIS = align_ephemeris(START_DATE, 52590., 0.09089812328, -1.59836e-13)

# Magnetar model parameters.
MAG_MODEL = xMagnetarModelsT2020.BLACKBODY
MAG_CHI = 85.
MAG_XI = 55.
MAG_DELTA_PHI = 0.35
MAG_BETA = 0.39

# Load the magnetar model table.
model_table = xMagnetarTableModelT2020()

# Interpolate the model table to the target parameters.
args = MAG_MODEL, MAG_CHI, MAG_XI, MAG_DELTA_PHI, MAG_BETA, INTEGRAL_FLUX
kwargs = dict(emin=INT_FLUX_EMIN, emax=INT_FLUX_EMAX)
spec, pol_deg, pol_ang = model_table.interpolate(*args, **kwargs)

# Define the actual ROI model.
psr = xPeriodicPointSource(SRC_NAME, SRC_RA, SRC_DEC, spec, pol_deg, pol_ang, EPHEMERIS, NH)
egb = xExtragalacticBkg(PNT_RA, PNT_DEC)
dge = xGalacticBkg(PNT_RA, PNT_DEC, GAL_BKG_RATE)
ROI_MODEL = xROIModel(PNT_RA, PNT_DEC, psr, egb, dge, instrumental_bkg)


def display_spectrum(emin=1., emax=12.):
    """Display the energy spectrum.
    """
    # Energy spectra at different phases.
    plt.figure('%s spectrum phase slices' % __model__)
    energy = numpy.linspace(emin, emax, 200)
    for phase in numpy.linspace(0., 1., 6):
        plt.plot(energy, spec(energy, phase), label='Phase = %.2f' % phase)
        plt.plot(energy, spec(6., phase) * (energy / 6.) ** -2.4, ls='dashed',
                 color=last_line_color())
    setup_gca(**fmtaxis.spec, xmin=emin, xmax=emax, logx=True, logy=True, grids=True, legend=True)
    # Pulse profiles at different energies.
    plt.figure('%s pulse profile' % __model__)
    phase = numpy.linspace(0., 1., 200)
    for energy in numpy.linspace(1., 8., 4):
        plt.plot(phase, spec(energy, phase), label='Energy = %.2f keV' % energy)
    setup_gca(**fmtaxis.pp_flux, grids=True, legend=True)


def display_pol_deg(emin=1., emax=12.):
    """Display the polarization degree.
    """
    # Slices of the polarization degree at different pulse-phase values.
    plt.figure('%s polarization degree phase slices' % __model__)
    energy = numpy.linspace(emin, emax, 200)
    for phase in numpy.linspace(0., 1., 6):
        plt.plot(energy, pol_deg(energy, phase), label='Phase = %.2f' % phase)
    setup_gca(**fmtaxis.ene_pol_deg, xmin=emin, xmax=emax, grids=True, legend=True)
    # Slices of the polarization degree at different energies.
    plt.figure('%s polarization degree pulse profile' % __model__)
    phase = numpy.linspace(0., 1., 200)
    for energy in numpy.linspace(1., 8., 4):
        plt.plot(phase, pol_deg(energy, phase), label='Energy = %.2f keV' % energy)
    setup_gca(**fmtaxis.pp_pol_deg, grids=True, legend=True)


def display_pol_ang(emin=1., emax=12.):
    """Display the polarization angle.
    """
    # Slices of the polarization angle at different pulse-phase values.
    plt.figure('%s polarization angle phase slices' % __model__)
    energy = numpy.linspace(emin, emax, 200)
    for phase in numpy.linspace(0., 1., 6):
        plt.plot(energy, numpy.degrees(pol_ang(energy, phase)), label='Phase = %.2f' % phase)
    setup_gca(**fmtaxis.ene_pol_ang, xmin=emin, xmax=emax, grids=True, legend=True)
    # Slices of the polarization degree at different energies.
    plt.figure('%s polarization angle pulse profile' % __model__)
    phase = numpy.linspace(0., 1., 200)
    for energy in numpy.linspace(1., 8., 4):
        plt.plot(phase, numpy.degrees(pol_ang(energy, phase)), label='Energy = %.2f keV' % energy)
    setup_gca(**fmtaxis.pp_pol_ang, grids=True, legend=True)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    display_spectrum(emin, emax)
    display_pol_deg(emin, emax)
    display_pol_ang(emin, emax)



if __name__ == '__main__':
    bootstrap_display()

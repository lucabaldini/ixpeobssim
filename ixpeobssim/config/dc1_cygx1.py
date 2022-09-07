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

Spectral Model
~~~~~~~~~~~~~~

The spectral model for Cygnus X-1 in the hard state is taken from
`Parker at al., 2015 <https://iopscience.iop.org/article/10.1088/0004-637X/808/1/9>`_
(Figure 6 and Table 4). In the IXPE energy range this includes three components:

* a comptonization component (eqpair);
* a reflection component (relxilllp);
* a narrow line at 6.4 keV.


.. _figure-cygx1_spectrum_e2:
.. figure:: figures/dc1/dc1_cygx1_spectrum_e2.*
   :width: 80%

.. _figure-cygx1_spectrum:
.. figure:: figures/dc1/dc1_cygx1_spectrum.*
   :width: 80%

   Input spectrum.



Polarization Model
~~~~~~~~~~~~~~~~~~

Giorgio suggested to use a 1--2% polarization degree for the main (eqpair) component,
with an arbitrary angle, and a 2% polarization degree, oriented at 90 degrees,
for the reflection component. The line should be unpolarized.

Also, he suggested not to include any temporal variability.


.. _figure-cygx1_polarization_degree:
.. figure:: figures/dc1/dc1_cygx1_polarization_degree.*
   :width: 80%

.. _figure-cygx1_polarization_angle:
.. figure:: figures/dc1/dc1_cygx1_polarization_angle.*
   :width: 80%

   Polarization model.


For completeness, here is an MDP table for the brightest source component.

.. code-block::

    xPointSource "eqpair" (id = 0)
        Galactic column density: 5.100e+21 cm^{-2}
        Redshift: 0.000
        Unabsorbed flux @ t = 0: 4.832e-09 erg/cm2/s (241.60 mcrab)
        Position: RA = 299.59031591 deg, Dec = 35.20160625 deg
    MDP table for 300000.0 s observation time
    2.00--2.83 keV:  4407218.0 counts, effective mu = 0.236, MDP = 0.87%
    2.83--4.00 keV:  2484315.3 counts, effective mu = 0.380, MDP = 0.72%
    4.00--5.66 keV:   891665.1 counts, effective mu = 0.473, MDP = 0.96%
    5.66--8.00 keV:   231725.5 counts, effective mu = 0.564, MDP = 1.58%
    2.00--8.00 keV:  8014923.9 counts, effective mu = 0.316, MDP = 0.48%

"""

import os

import numpy

from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII
from ixpeobssim.config import file_path_to_model_name, bootstrap_display
from ixpeobssim.config.dc1_bkg import instrumental_bkg
from ixpeobssim.config.dc1_utils import sum_point_sources
from ixpeobssim.core.spline import xInterpolatedUnivariateLogSpline
from ixpeobssim.srcmodel.bkg import xExtragalacticBkg, xGalacticBkg
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import gaussian_line
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.units_ import ergcms_to_mcrab


#pylint: disable=invalid-name


__model__ = file_path_to_model_name(__file__)


# Source coordinates, in decimal degrees.
SRC_NAME = 'Cyg X-1'
SRC_RA, SRC_DEC = 299.59031591, 35.20160625
SRC_L, SRC_B = 71.33499793, 3.0668341

# Pointing coordinates
PNT_RA, PNT_DEC = SRC_RA, SRC_DEC

# Observation start date
START_DATE = '2022-04-01'
DURATION = 550000.

# Interstellar absorption.
NH = 5.1e21

# Galactic diffuse.
GAL_BKG_RATE = 70.

# Polarization patterns for the various components.
EQPAIR_POL_DEG = 0.025
EQPAIR_POL_ANG = 30.

RELXILLLP_POL_DEG = 0.020
RELXILLLP_POL_ANG = EQPAIR_POL_ANG - 90.

LINE_POL_DEG = 0.
LINE_POL_ANG = 0.


def _parse_spec(file_name):
    """Parse the spectral data.

    In absence of the component normalization in Parker et al., 2015, we extracted
    this data from Figure 6 (right panel) of the paper, to be used as a reference.

    Note that only the eqpair and relxilllp components have been digitized.
    """
    file_path = os.path.join(IXPEOBSSIM_CONFIG_ASCII, file_name)
    logger.info('Loading spectral data from %s...', file_path)
    energy, e2flux = numpy.loadtxt(file_path, unpack=True, delimiter=',')
    flux = e2flux / (energy**2.)
    spec_spline = xInterpolatedUnivariateLogSpline(energy, flux, k=1, **fmtaxis.spec)
    e2spec_spline = xInterpolatedUnivariateLogSpline(energy, e2flux, k=1, **fmtaxis.spec)
    return spec_spline, e2spec_spline


eqpair_spec_spline, eqpair_e2spec_spline = _parse_spec('cyg_x1_eqpair_e2spec.csv')
eqpair_spec = lambda E, t=None: eqpair_spec_spline(E)
eqpair_pol_deg = constant(EQPAIR_POL_DEG)
eqpair_pol_ang = constant(numpy.radians(EQPAIR_POL_ANG))
eqpair = xPointSource('eqpair', SRC_RA, SRC_DEC, eqpair_spec, eqpair_pol_deg,
                      eqpair_pol_ang, column_density=NH)

relxilllp_spec_spline, relxilllp_e2spec_spline = _parse_spec('cyg_x1_relxilllp_e2spec.csv')
relxilllp_spec = lambda E, t=None: relxilllp_spec_spline(E)
relxilllp_pol_deg = constant(RELXILLLP_POL_DEG)
relxilllp_pol_ang = constant(numpy.radians(RELXILLLP_POL_ANG))
relxilllp = xPointSource('relxilllp', SRC_RA, SRC_DEC, relxilllp_spec,
                         relxilllp_pol_deg, relxilllp_pol_ang, column_density=NH)

line_spec = gaussian_line(1.e-3, 6.4, 0.01)
line_pol_deg = constant(LINE_POL_DEG)
line_pol_ang = constant(numpy.radians(LINE_POL_ANG))
line = xPointSource('line', SRC_RA, SRC_DEC, line_spec, line_pol_deg, line_pol_ang,
                    column_density=NH)

e, total_spec, total_pol_deg, total_pol_ang = sum_point_sources(eqpair, relxilllp, line)

# Create the background components.
egb = xExtragalacticBkg(PNT_RA, PNT_DEC)
dge = xGalacticBkg(PNT_RA, PNT_DEC, GAL_BKG_RATE)

# Create the ROI model.
ROI_MODEL = xROIModel(PNT_RA, PNT_DEC, eqpair, relxilllp, line, egb, dge, instrumental_bkg)


def display_spectrum(energy):
    """Plot the energy spectrum of the source.
    """
    # Calculate the standard integral flux.
    emin = 2.
    emax = 10.
    source_dict = {name: src for name, src in ROI_MODEL.items() if src.identifier < 3}
    flux_ergcms = sum(src.calculate_integral_flux(emin, emax, t=0) for src in source_dict.values())
    flux_mcrab = ergcms_to_mcrab(flux_ergcms)
    flux_label = 'Total flux: %.2e erg cm$^{-2}$ s$^{-1}$ (%.1f mcrab)' % (flux_ergcms, flux_mcrab)

    plt.figure('%s spectrum E2' % __model__)
    for name, component in source_dict.items():
        plt.plot(energy, energy**2. * component.photon_spectrum(energy), label=name)
    plt.plot(e, e**2. * total_spec, label='Total')
    plt.text(1.2, 7.5, flux_label)
    setup_gca(xmin=energy.min(), xmax=energy.max(), ymin=1.e-1, logx=True,
              logy=True, grids=True, legend=True, **fmtaxis.spec)

    plt.figure('%s spectrum' % __model__)
    for name, component in source_dict.items():
        plt.plot(energy, component.photon_spectrum(energy), label=name)
    plt.plot(e, total_spec, label='Total')
    setup_gca(xmin=energy.min(), xmax=energy.max(), ymin=1.e-3, logx=True,
              logy=True, grids=True, legend=True, **fmtaxis.spec)


def display_polarization(energy):
    """Plot the polarization degree and angle as a function of the energy.
    """
    source_dict = {name: src for name, src in ROI_MODEL.items() if src.identifier < 3}
    plt.figure('%s polarization degree' % __model__)
    for name, component in source_dict.items():
        plt.plot(energy, component.polarization_degree(energy), label=name)
    plt.plot(e, total_pol_deg, label='Total')
    setup_gca(xmin=energy.min(), xmax=energy.max(), ymin=-0.001, ymax=0.03,
              legend=True, grids=True, **fmtaxis.ene_pol_deg)

    plt.figure('%s polarization angle' % __model__)
    for name, component in source_dict.items():
        plt.plot(energy, numpy.degrees(component.polarization_angle(energy)), label=name)
    plt.plot(e, numpy.degrees(total_pol_ang), label='Total')
    setup_gca(xmin=energy.min(), xmax=energy.max(), ymin=-90., ymax=90.,
              legend=True, grids=True, **fmtaxis.ene_pol_ang_deg)


def display(emin=1., emax=15.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 1000)
    display_spectrum(energy)
    display_polarization(energy)



if __name__ == '__main__':
    bootstrap_display()

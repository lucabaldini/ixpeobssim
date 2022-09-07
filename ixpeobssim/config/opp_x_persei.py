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

"""Simple model for X - Persei.

Ephemeris taken from H. Delgado-Martí et al. (2001)

Pulse profile and spectrum taken from Di Salvo et al. (2009). For the pulse
profile we use fig. 3, middle panel, corresponding to the energy range
bewteen 1.8 and 10.5 keV.
"""

import os

import numpy

from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII
from ixpeobssim.config import file_path_to_model_name, bootstrap_display
from ixpeobssim.core.spline import xUnivariateSpline
from ixpeobssim.srcmodel.ephemeris import xEphemeris
from ixpeobssim.srcmodel.roi import xPeriodicPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import highecut_power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.utils.matplotlib_ import plt, setup_gca


__model__ = file_path_to_model_name(__file__)


# Source coordinates, in decimal degrees.
SRC_NAME = 'X Persei'
SRC_RA, SRC_DEC = 58.84615783, 31.04584604
SRC_L, SRC_B = 163.08135735, -17.13618855

# Pointing coordinates
PNT_RA, PNT_DEC = SRC_RA, SRC_DEC

# Interstellar absorption.
NH = 0.151e22

# Source spectrum.
HECPL_NORM = 0.015
HECPL_INDEX = 0.441
HECPL_ECUT = 2.180
HECPL_EFOLD = 3.94

# Source ephemeris.
EPHEMERIS = xEphemeris(0., 837.6713)

# Source polarization
POL_DEG = 0.05
POL_ANG = 0.


def _parse_pulse_profile():
    """Parse the X Persei spectral data from a csv file.

    The data were extracted from Figure 3, mid panel from Di Salvo et al. (1998).

    Note we don't use the phase value in the file, as the underlying data points
    have 64 equi-spaced bins, and we can calculate the bin centers without
    relying on the manual image processing.
    """
    file_path = os.path.join(IXPEOBSSIM_CONFIG_ASCII, 'xpersei_pp.csv')
    x, y = numpy.loadtxt(file_path, delimiter=',', unpack=True)
    num_bins = 64
    x = numpy.linspace(0., 1. - 1. / num_bins, num_bins) + 0.5 / num_bins
    # Normalize the y values so that the average is one.
    y /= y.mean()
    return x, y


pp_phase, pp_norm = _parse_pulse_profile()
pulse_profile = xUnivariateSpline(pp_phase, pp_norm, s=0.025)
spec_norm = lambda t: HECPL_NORM * pulse_profile(t)
spec = highecut_power_law(spec_norm, HECPL_INDEX, HECPL_ECUT, HECPL_EFOLD)
pol_deg = constant(POL_DEG)
pol_ang = constant(POL_ANG)
src = xPeriodicPointSource(SRC_NAME, SRC_RA, SRC_DEC, spec, pol_deg, pol_ang,
                          ephemeris=EPHEMERIS, column_density=NH)

ROI_MODEL = xROIModel(PNT_RA, PNT_DEC, src)


def display(emin=1., emax=12.):
    """Main display() function.
    """
    energy = numpy.linspace(emin, emax, 200)
    plt.figure('%s spectrum' % __model__)
    for phase in numpy.arange(0., 1., 0.2):
        plt.plot(energy, spec(energy, phase), label='Phase = %.2f' % phase)
    setup_gca(xmin=energy.min(), xmax=energy.max(), logx=True, logy=True,
              grids=True, legend=True, **fmtaxis.spec)

    plt.figure('%s pulse profile' % __model__)
    plt.plot(pp_phase, pp_norm, 'o', label='Input data')
    pulse_profile.plot(label='Spline')
    setup_gca(grids=True, legend=True, xlabel='Pulse phase',
              ylabel='Spectrum normalization')



if __name__ == '__main__':
    bootstrap_display()

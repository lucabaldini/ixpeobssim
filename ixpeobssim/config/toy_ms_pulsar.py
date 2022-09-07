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

"""Toy ms pulsar with a phase-dependent energy spectrum.

This was specifically design to study the deadtime correction, the main things
in mind being:

* we want a period comparable or smaller than the average dead time per event,
to make the correction maximally difficult;
* we want the normalization of the energy spectrum to span all the different
regimes (from faint to bright) through the pulse phase;
* we want a significant spectral evolution throuough the phase so that the
average energy changes measurably---and so does the average ROI size, which
in turns determines the dead time per event.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.srcmodel.roi import xPeriodicPointSource, xROIModel
from ixpeobssim.srcmodel.ephemeris import xEphemeris
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.utils.matplotlib_ import plt, setup_gca


__model__ = file_path_to_model_name(__file__)

SRC_RA, SRC_DEC = 45., 30.
BASELINE_NORM = 0.2
BASELINE_INDEX = 3.
PEAK1_NORM_MAX = 12.
PEAK1_CENTER = 0.2
PEAK1_WIDTH = 0.055
PEAK1_INDEX_HARDENING = -1.5
PEAK2_NORM_MAX = 4.
PEAK2_CENTER = 0.8
PEAK2_WIDTH = 0.060
PEAK2_INDEX_HARDENING = -1.
POL_DEG = 0.3
POL_ANG = 40.
NU_0 = 200.


def _gaussian_peak(phase, max, center, width):
    """Simple Gaussian paramtrization for the purpose of generating multiple-peak
    synthetic pulse profiles.
    """
    return max * numpy.exp(-(phase - center)**2. / (2. * width**2.))

def pl_norm(phase):
    """Power-law normalization as a function of the pulse phase.
    """
    return BASELINE_NORM + \
        _gaussian_peak(phase, PEAK1_NORM_MAX, PEAK1_CENTER, PEAK1_WIDTH) +\
        _gaussian_peak(phase, PEAK2_NORM_MAX, PEAK2_CENTER, PEAK2_WIDTH)

def pl_index(phase):
    """Power-law index as a function of the pulse phase.
    """
    return BASELINE_INDEX + \
        _gaussian_peak(phase, PEAK1_INDEX_HARDENING, PEAK1_CENTER, PEAK1_WIDTH) +\
        _gaussian_peak(phase, PEAK2_INDEX_HARDENING, PEAK2_CENTER, PEAK2_WIDTH)


spec = power_law(pl_norm, pl_index)
pol_deg = constant(POL_DEG)
pol_ang = constant(numpy.radians(POL_ANG))

ephemeris = xEphemeris(0, NU_0)
src = xPeriodicPointSource('Periodic source', SRC_RA, SRC_DEC, spec, pol_deg, pol_ang, ephemeris)
src.set_count_spectrum_params(500, 3, 1)
ROI_MODEL = xROIModel(SRC_RA, SRC_DEC, src)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    phase = numpy.linspace(0., 1., 200)
    energy = numpy.linspace(emin, emax, 200)
    plt.figure('%s PL normalization' % __model__)
    plt.plot(phase, pl_norm(phase))
    setup_gca(grids=True, **fmtaxis.pp_pl_norm)
    plt.figure('%s PL index' % __model__)
    plt.plot(phase, pl_index(phase))
    setup_gca(grids=True, ymax=BASELINE_INDEX + 0.5, **fmtaxis.pp_pl_index)
    plt.figure('%s spectral shape' % __model__)
    for p in (PEAK1_CENTER, 0.5, PEAK2_CENTER):
        plt.plot(energy, spec(energy, p), label='Phase = %.3f' % p)
    setup_gca(xmin=emin, xmax=emax, logx=True, logy=True, legend=True, grids=True)
    plt.show()



if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

# Copyright (C) 2022, the ixpeobssim team.
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

"""Minimal example created for the SoftwareX paper.

The field contains a single, stationary point source with a power-law photon
spectrum and a polarization degree linearly increasing with energy (with the
corresponding position angle constant, and aligned with the celestial North).
"""

import numpy

from ixpeobssim.srcmodel.bkg import xTemplateInstrumentalBkg
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law

# Sky position and spectral parameters of the source.
SRC_RA, SRC_DEC = 20., 30.
PL_NORM = 6.e-3
PL_INDEX = 2.
# The pointing direction is the same as the source coordinates.
PNT_RA, PNT_DEC = SRC_RA, SRC_DEC

# Definition of the photon spectrum.
spec = power_law(PL_NORM, PL_INDEX)

def pol_deg(E, t=None, ra=None, dec=None):
    """Definition of the polarization degree as a function of the energy.

    The polarization degree is 5% at 1 keV, increasing linearly with energy.
    """
    return 0.05 * E

# Definition of the polarization angle---0. is aligned with the North.
pol_ang = constant(0.)

# Definition of the sources and the region of interest.
src = xPointSource('Point source', SRC_RA, SRC_DEC, spec, pol_deg, pol_ang)
bkg = xTemplateInstrumentalBkg()
ROI_MODEL = xROIModel(PNT_RA, PNT_DEC, src, bkg)


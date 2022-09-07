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

"""Toy point source with (non-negligible) instrumental background.

This is useful to study the background subtraction in spectro-polarimetric analysis.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name, bootstrap_display
from ixpeobssim.srcmodel.bkg import xTemplateInstrumentalBkg
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.utils.matplotlib_ import plt, setup_gca


__model__ = file_path_to_model_name(__file__)


# Source coordinates, in decimal degrees.
SRC_NAME = 'Toy point surce w/ bkg'
SRC_RA, SRC_DEC = 45., 45.

# Pointing coordinates
PNT_RA, PNT_DEC = SRC_RA, SRC_DEC

# Spectral and polarimetric parameters
PL_NORM = 6.e-3
PL_INDEX = 3.
PD = 0.5
PA = 30.
SPEC = power_law(PL_NORM, PL_INDEX)
POL_DEG = constant(PD)
POL_ANG = constant(numpy.radians(PA))

# Definition of the sources and the region of interest.
SRC = xPointSource('Point source', SRC_RA, SRC_DEC, SPEC, POL_DEG, POL_ANG)
BKG = xTemplateInstrumentalBkg()
ROI_MODEL = xROIModel(PNT_RA, PNT_DEC, SRC, BKG)


def display(emin=1., emax=12., area=3.e-2):
    """Display the source model.

    Note the background needs to be multiplied for an effective area on the
    readout plane to be compared with the source. A circle with 1~mm radius is
    a sensible proxy for a source cut.
    """
    energy = numpy.linspace(emin, emax, 100)
    plt.figure('%s spectrum' % __model__)
    plt.plot(energy, SRC.photon_spectrum(energy), label='Source')
    plt.plot(energy, area * BKG.photon_spectrum(energy), label='Background')
    setup_gca(xmin=emin, xmax=emax, logx=True, logy=True, grids=True, legend=True, **fmtaxis.spec)


if __name__ == '__main__':
    bootstrap_display()

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

"""Simple model for Mrk 501.

The spectral model is adapted from B. Kapanadze et al., MNRAS469,1655-1672 (2017)
"""

import numpy

from ixpeobssim.config import file_path_to_model_name, bootstrap_display
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.utils.matplotlib_ import plt, setup_gca


__model__ = file_path_to_model_name(__file__)


# Source coordinates, in decimal degrees.
SRC_NAME = 'Mrk 501'
SRC_RA, SRC_DEC = 253.46756952, 39.76016915
SRC_L, SRC_B = 63.60003969, 38.85915565

# Pointing coordinates
PNT_RA, PNT_DEC = SRC_RA, SRC_DEC

# Interstellar absorption.
NH = 1.55e20

# Source spectrum.
PL_NORM = 0.05
PL_INDEX = 1.9

# Source polarization
POL_DEG = 0.05
POL_ANG = 0.

spec = power_law(PL_NORM, PL_INDEX)
pol_deg = constant(POL_DEG)
pol_ang = constant(POL_ANG)
src = xPointSource(SRC_NAME, SRC_RA, SRC_DEC, spec, pol_deg, pol_ang, column_density=NH)

ROI_MODEL = xROIModel(PNT_RA, PNT_DEC, src)


def display(emin=1., emax=12.):
    """Main display() function.
    """
    energy = numpy.linspace(emin, emax, 200)
    plt.figure('%s spectrum' % __model__)
    plt.plot(energy, spec(energy))
    setup_gca(xmin=energy.min(), xmax=energy.max(), logx=True, logy=True,
              grids=True, legend=True, **fmtaxis.spec)



if __name__ == '__main__':
    bootstrap_display()

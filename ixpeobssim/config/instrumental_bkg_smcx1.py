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

"""First instrumental background parametrization based on actual data---a short
observation of SMC X-1.
"""

from __future__ import print_function, division

import os

import numpy

from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xROIModel
from ixpeobssim.srcmodel.bkg import xTemplateInstrumentalBkg


__model__ = file_path_to_model_name(__file__)

RA = 45.
DEC = 45.

bkg = xTemplateInstrumentalBkg()

ROI_MODEL = xROIModel(RA, DEC, bkg)


def display(emin=0., emax=15.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 250)
    # Energy spectrum
    plt.figure('%s spectrum' % __model__)
    plt.plot(energy, bkg.photon_spectrum(energy))
    setup_gca(xmin=emin, xmax=emax, ymin=bkg.photon_spectrum(emax),
              logy=True, grids=True, **fmtaxis.spec)



if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

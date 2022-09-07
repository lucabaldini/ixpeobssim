#!/usr/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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

"""Sample configuration file for the celestial background.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xROIModel
from ixpeobssim.srcmodel.bkg import xExtragalacticBkg, xGalacticBkg


__model__ = file_path_to_model_name(__file__)

ra = 45.
dec = 45.


egb = xExtragalacticBkg(ra, dec)
dge = xGalacticBkg(ra, dec, 1000.)

ROI_MODEL = xROIModel(ra, dec, egb, dge)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 100)
    # Energy spectrum
    plt.figure('%s spectrum' % __model__)
    plt.plot(energy, egb.photon_spectrum(energy), label='Extragalactic background')
    plt.plot(energy, dge.photon_spectrum(energy), label='Galactic background')
    setup_gca(xmin=emin, xmax=emax, ymin=egb.photon_spectrum(emax), legend=True,
              logx=True, logy=True, grids=True, **fmtaxis.spec)



if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

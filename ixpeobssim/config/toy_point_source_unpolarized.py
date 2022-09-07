#!/usr/bin/env python
#
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

"""Unpolarized point source with a power-law spectrum.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.fmtaxis import fmtaxis


__model__ = file_path_to_model_name(__file__)
ra, dec = 30., 45.
pl_norm = 1.
pl_index = 2.
spec = power_law(pl_norm, pl_index)
pol_deg = constant(0.)
pol_ang = constant(0.)

src = xPointSource('Point source', ra, dec, spec, pol_deg, pol_ang)

ROI_MODEL = xROIModel(ra, dec, src)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 100)

    # Energy spectrum
    plt.figure('%s spectrum' % __model__)
    plt.plot(energy, spec(energy))
    setup_gca(xmin=emin, xmax=emax, ymin=spec(emax), logx=True, logy=True,
              grids=True, **fmtaxis.spec)


if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

#!/usr/bin/env python
#
# Copyright (C) 2020, the ixpeobssim team.
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

"""
This is a point source with a time-independent power-law spectrum. The
polarization degrees is linearly increasing with energy, and the polarization
angle is constant.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.fmtaxis import fmtaxis

# pylint: disable=invalid-name

__model__ = file_path_to_model_name(__file__)
ra, dec = 45., 45.
pl_norm = 10.
pl_index = 2.
pd_intercept = 0.1
pd_slope = 0.1
spec = power_law(pl_norm, pl_index)

def pol_deg(E, t=None, ra=None, dec=None):
    """Linear model for the polarization degree vs. energy.
    """
    return numpy.clip(pd_intercept + pd_slope * (E - 1.), 0., 1.)

pa = 30.
pol_ang = constant(numpy.radians(pa))

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
    # Polarization degree.
    plt.figure('%s polarization degree' % __model__)
    plt.plot(energy, pol_deg(energy))
    setup_gca(xmin=emin, xmax=emax, ymin=0, grids=True, **fmtaxis.ene_pol_deg)



if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

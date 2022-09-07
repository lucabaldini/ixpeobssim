#!/usr/bin/env python
#
# Copyright (C) 2015--2020, the ixpeobssim team.
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
This is a simple example of a point source with a complex spectral model from
XSPEC. (Note this can be used along with any source morphology or polarization
pattern.)

The only meaningful feature of the source that can be plotted is its
energy spectrum, which is shown below.

.. image:: figures/models/toy_xspec_point_source_spectrum.png
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.srcmodel.spectrum import xXspecModel


__model__ = file_path_to_model_name(__file__)
ra, dec = 45., 45.
# Parameters for a complex XSPEC spectral model phabs*(bbody+powerlaw)
expression = 'phabs * (bbody + powerlaw)'
col_dens = 1.
bb_temp = 1. #in keV
# bb_norm is the luminosity in units of 1e39 erg/s over
# the distance squared in units of 10 kpc
bb_norm = 1e-3 / (0.1**2)
pl_index = 2.
pl_norm = 10.

spec = xXspecModel(expression, [col_dens, bb_temp, bb_norm, pl_index, pl_norm])
pd = 0.5
pa = 30.
pol_deg = constant(pd)
pol_ang = constant(numpy.radians(pa))

src = xPointSource('Point source', ra, dec, spec, pol_deg, pol_ang)

ROI_MODEL = xROIModel(ra, dec, src)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    plt.figure('%s spectrum' % __model__)
    spec.plot()
    setup_gca(xmin=emin, xmax=emax, ymin=0.1, ymax=3., logx=True, logy=True, grids=True)



if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

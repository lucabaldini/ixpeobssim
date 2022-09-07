#!/usr/bin/env python
#
# Copyright (C) 2015--2018, the ixpeobssim team.
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
This is very similar to the uniform disk example, except that the morpholgy
of the source is a bivariate gaussian (with a diagonal covariance matrix)
and the source does not extend as much across the field of view, so that it 
is pretty much completely contained.

Below is yet another binned coount map from xpbin in CMAP mode.

.. image:: figures/obssim/toy_gauss_disk_count_map.png
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xGaussianDisk, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.units_ import arcmin_to_degrees
from ixpeobssim.utils.fmtaxis import fmtaxis


__model__ = file_path_to_model_name(__file__)
ra, dec = 45., 45.
radius = arcmin_to_degrees(2.) 
pl_norm = 1.
pl_index = 2.
spec = power_law(pl_norm, pl_index)
pol_deg = constant(0.5)
pol_ang = constant(numpy.radians(65.))

src = xGaussianDisk('Gaussian disk', ra, dec, radius, spec, pol_deg, pol_ang)

ROI_MODEL = xROIModel(ra, dec, src)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 100)

    # Energy spectrum
    plt.figure('%s spectrum' % __model__)
    plt.plot(energy, spec(energy), label=src.name)
    setup_gca(xmin=emin, xmax=emax, ymin=spec(emax), logx=True, logy=True,
              grids=True, **fmtaxis.spec)


if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

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
This is possibly the simplest interesting model that can be simulated: a
point source with a time-independent power-law spectrum and time- and
energy-independent polarization degree and angle.

The only meaningful feature of the source that can be plotted is its
energy spectrum, which is shown below.

.. image:: figures/models/toy_point_source_spectrum.png

And here is a spectral fit with XSPEC to a set of event lists created
with xpobssim and binned with xpbin:

.. image:: figures/obssim/toy_point_source_xspec_fit.png
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
pl_norm = 10.
pl_index = 2.
spec = power_law(pl_norm, pl_index)
pd = 0.1
pa = 30.
pol_deg = constant(pd)
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


if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

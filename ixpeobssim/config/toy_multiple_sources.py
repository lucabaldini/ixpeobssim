#!/usr/bin/env python
#
# Copyright (C) 2015--2016, the ixpeobssim team.
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
This is a slight variation over the single point source, where we include
in the same region of interest two stationary point sources, with different
power-law spectra and polarization patterns, offset by a few arcminutes with
respect to the center of the ROI.

For completeness, below are the energy spectra for the two sources.

.. image:: figures/models/toy_multiple_sources_spectrum.png

And here is a count map created from a simulated observation of this ROI model
using xpbin in CMAP mode.

.. image:: figures/obssim/toy_multiple_sources_count_map.png

.. tip::
   In this case we are using the ability of ixpeobssim to overlay multiple
   sources in the same ROI to simulate two physically distinct objects, but
   subclasses of :class:`ixpeobssim.srcmodel.roi.xModelComponentBase` should
   in fact be thought of as `model components`, and this is the very same
   technique that we use to simulate different physical components of the same
   object.
"""

from __future__ import print_function, division


import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.units_ import arcmin_to_degrees
from ixpeobssim.utils.fmtaxis import fmtaxis


__model__ = file_path_to_model_name(__file__)
ra, dec = 30., 45.

ra1, dec1 = ra, dec
pl_norm1 = 1.
pl_index1 = 2.
spec1 = power_law(pl_norm1, pl_index1)
pol_deg1 = constant(0.25)
pol_ang1 = constant(numpy.radians(30.))

ra2, dec2 = ra + arcmin_to_degrees(2.), dec - arcmin_to_degrees(3.)
pl_norm2 = 0.25
pl_index2 = 1.5
spec2 = power_law(pl_norm2, pl_index2)
pol_deg2 = constant(0.)
pol_ang2 = constant(numpy.radians(0.))

src1 = xPointSource('Point source 1', ra1, dec1, spec1, pol_deg1, pol_ang1)
src2 = xPointSource('Point source 2', ra2, dec2, spec2, pol_deg2, pol_ang2)

ROI_MODEL = xROIModel(ra, dec, src1, src2)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 100)

    # Energy spectrum
    plt.figure('%s spectrum' % __model__)
    plt.plot(energy, spec1(energy), label=src1.name)
    plt.plot(energy, spec2(energy), label=src2.name)
    setup_gca(xmin=emin, xmax=emax, ymin=spec2(emax), logx=True, logy=True,
              grids=True, legend=True, **fmtaxis.spec)


if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

#!/usr/bin/env python
#
# Copyright (C) 2018, the ixpeobssim team.
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
This is a fairly complicated toy example, illustrating a number of ixpeobssim
features.

The ROI includes a single point source with a power-law spectrum, but this
time around, the source being periodic, the corresponding power-law
normalization varies with the pulse phase according to a simple analytical
formula. (In order to make the comparison of the simulation output with the
input model straightforward the power-law index is constant.) Below is the
normalization as a function of the pulse phase

.. image:: figures/models/toy_periodic_source_pl_normalization.png

along with the actual differential energy spectrum, plotted at several
different pulse phases.

.. image:: figures/models/toy_periodic_source_spectrum.png

The polarization angle is constant, but for illustration purposes the
polarization degree varies with the pulse phase according to yet another
simple formula, as shown in the plot below:

.. image:: figures/models/toy_periodic_source_polarization_degree.png

Here is an example of a binned pulse profile from an event list
produced by ixpeobssim, compared with the input model.

.. image:: figures/obssim/toy_periodic_source_pulse_profile.png

And below is an example of a slightly more complicated analysis, where we

* fold the original event list with the ephemeris of the input model;
* split the original even list into a number of phase bins;
* estimate the polarization degree into two energy bins (within each phase bin)
  and compare with the input model.

.. image:: figures/obssim/toy_periodic_source_polarization_degree.png

For completeness, the polarization degree is evaluated by summing the
event-by-event Stokes parameters, i.e., without taking into account the
energy redstribution---which is probably account for most of the slight bias
that can be see in the plot. The input model is evaluated at the average
photon energy within each bin.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xPeriodicPointSource, xROIModel
from ixpeobssim.srcmodel.ephemeris import xEphemeris
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.utils.time_ import string_to_met_utc


__model__ = file_path_to_model_name(__file__)
ra, dec = 45., 45.

def pl_norm(phase):
    """Energy spectrum: power-law normalization as a function of the pulse
    phase.
    """
    return 2. * (1.25 + numpy.cos(4 * numpy.pi * phase))

pl_index = 2.
spec = power_law(pl_norm, pl_index)

def pol_deg(E, phase, ra=None, dec=None):
    """Polarization degree as a function of the dynamical variables.

    Since we're dealing with a point source the sky direction (ra, dec) is
    irrelevant and, as they are not used, defaulting the corresponding arguments
    to None allows to call the function passing the energy and phase only.
    """
    norm = numpy.clip(E / 10., 0., 1.)
    return norm * (0.5 + 0.25 * numpy.cos(4 * numpy.pi * (phase - 0.25)))

pol_ang = constant(numpy.radians(65.))

start_date = '2022-04-21'
t0 = string_to_met_utc(start_date, lazy=True)
nu0 = 0.1363729749
nudot0 = 0.543638e-13
ephemeris = xEphemeris(t0, nu0, nudot0)
src = xPeriodicPointSource('Periodic source', ra, dec, spec, pol_deg, pol_ang, ephemeris)

ROI_MODEL = xROIModel(ra, dec, src)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 100)
    phase = numpy.linspace(0., 1., 100)

    # Pulse profile: power-law normalization.
    plt.figure('%s PL normalization' % __model__)
    plt.plot(phase, pl_norm(phase))
    setup_gca(ymin=0., ymax=2.5, **fmtaxis.pp_pl_norm)

    # Pulse profile: polarization degree.
    plt.figure('%s polarization degree' % __model__)
    for E in [2., 5., 8.]:
        plt.plot(phase, pol_deg(E, phase), label='Energy = %.2f keV' % E)
    setup_gca(ymin=0., ymax=1., legend=True, **fmtaxis.pp_pol_deg)

    # Energy spectrum at different phase values.
    plt.figure('%s spectrum' % __model__)
    for p in [0.00, 0.05, 0.10, 0.15, 0.20, 0.25]:
        plt.plot(energy, spec(energy, p), label='Phase = %.2f' % p)
    setup_gca(xmin=emin, xmax=emax, logx=True, logy=True, legend=True,
              grids=True, **fmtaxis.spec)


if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

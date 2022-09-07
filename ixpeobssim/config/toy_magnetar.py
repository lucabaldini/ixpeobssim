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
This is an ad-hoc source to debug
https://bitbucket.org/ixpesw/ixpeobssim/issues/265/
(aka the magnetar issue).
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


__model__ = file_path_to_model_name(__file__)
source_name = 'Toy magnetar'
ra, dec = 257.2042, -40.1528
nu0 = 0.09089812328
pl_norm = 1.
pl_index = 2.
spec = power_law(pl_norm, pl_index)


def pol_deg(E, phase, ra=None, dec=None):
    """Polarization degree as a function of the dynamical variables.
    """
    return 0.5 + 0.25 * numpy.cos(4 * numpy.pi * (phase - 0.25))


def pol_ang(E, phase, ra=None, dec=None):
    """Polarization angle as a function of the dynamical variables.
    """
    return numpy.radians(30 + 18. * numpy.cos(4 * numpy.pi * (phase - 0.25)))


ROI_MODEL = xROIModel(ra, dec)
ephemeris = xEphemeris(0., nu0)
src = xPeriodicPointSource(source_name, ra, dec, spec, pol_deg, pol_ang, ephemeris)
ROI_MODEL.add_source(src)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 100)
    phase = numpy.linspace(0., 1., 100)

    # Energy spectrum.
    plt.figure('%s spectrum' % __model__)
    plt.plot(energy, spec(energy))
    setup_gca(xmin=emin, xmax=emax, logx=True, logy=True, grids=True, **fmtaxis.spec)

    # Pulse profile: polarization degree.
    plt.figure('%s polarization degree' % __model__)
    plt.plot(phase, pol_deg(None, phase))
    setup_gca(ymin=0., ymax=1., grids=True, **fmtaxis.pp_pol_deg)

    # Pulse profile: polarization angle.
    plt.figure('%s polarization angle' % __model__)
    plt.plot(phase, numpy.degrees(pol_ang(None, phase)))
    setup_gca(ymin=0., ymax=50., grids=True, **fmtaxis.pp_pol_ang)



if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

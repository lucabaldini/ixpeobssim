# Copyright (C) 2021, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public Licensese as published by
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


"""Temporary utils for the data challenge 1 models.
"""

import numpy

from ixpeobssim.core.stokes import xModelStokesParameters
from ixpeobssim.srcmodel.ephemeris import xEphemeris
from ixpeobssim.utils.time_ import mjd_to_met, string_to_met_utc


#pylint: disable=invalid-name


def sum_pol_deg_ang(*items):
    """Small throw-away convenience function.
    """
    I = sum(f for f, _, _ in items)
    Q = sum(f * xModelStokesParameters.q(pd, pa) for f, pd, pa in items)
    U = sum(f * xModelStokesParameters.u(pd, pa) for f, pd, pa in items)
    q = xModelStokesParameters.normalize(Q, I)
    u = xModelStokesParameters.normalize(U, I)
    pol_deg = xModelStokesParameters.polarization_degree(q, u)
    pol_ang = xModelStokesParameters.polarization_angle(q, u)
    return I, pol_deg, pol_ang


def sum_point_sources(*sources):
    """Quick and dirty point-source addition.

    This should be included in ixpeobssim, see issue #411
    """
    energy = numpy.linspace(1., 15., 1000)
    I = sum(src.photon_spectrum(energy) for src in sources)
    Q = sum(src.photon_spectrum(energy) * \
        xModelStokesParameters.q(src.polarization_degree(energy),\
        src.polarization_angle(energy)) for src in sources)
    U = sum(src.photon_spectrum(energy) * \
        xModelStokesParameters.u(src.polarization_degree(energy),\
        src.polarization_angle(energy)) for src in sources)
    q = xModelStokesParameters.normalize(Q, I)
    u = xModelStokesParameters.normalize(U, I)
    pol_deg = xModelStokesParameters.polarization_degree(q, u)
    pol_ang = xModelStokesParameters.polarization_angle(q, u)
    return energy, I, pol_deg, pol_ang


def align_ephemeris(start_date, t0, nu0, nudot0, nuddot=0):
    """Align a set of ephemeris to a given start date.

    Note this would be a conveniente addition to the xEphemeris class as a
    method.
    """
    start_met = string_to_met_utc(start_date, lazy=True)
    met0 = mjd_to_met(t0)
    dt = start_met - met0
    nu = nu0 + nudot0 * dt + nuddot * (dt**2.) / 2.
    nudot = nudot0 + nuddot * dt
    return xEphemeris(start_met, nu, nudot, nuddot)

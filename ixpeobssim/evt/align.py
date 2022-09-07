#!/usr/bin/env python
#
# Copyright (C) 2021 the ixpeobssim team.
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

"""Facilities for the rotation of the polarization angle/Stokes parameters
given an input model.

This module provides a series of functions facilitating the search for large-scale
polarization signatures, such as radial/tangential polarization in extended
sources such SNRs.
"""

from __future__ import print_function, division

from ixpeobssim.evt.kislat2015 import xStokesAnalysis
from ixpeobssim.utils.math_ import modulo_2pi


# pylint: disable=invalid-name


def align_phi(phi, phi0):
    """Rotate the photoelectron azimuthal angle and recalculate the event-by-event
    Stokes parameters.

    Args
    ----
    phi : array_like
        The azimuthal angle for the input event list

    phi0 : array_like
        The model polarization direction, calculated at the positions of the
        input events.

    .. warning::

        As we moved away from the used of azimuthal angles in the analysis, this
        function was added to the module for backward compatibility and testing
        purposes, but should no be used. Use align_stokes_parameters() below,
        instead.
    """
    aligned_phi = modulo_2pi(phi - phi0)
    aligned_q = xStokesAnalysis.stokes_q(aligned_phi, weights=None)
    aligned_u = xStokesAnalysis.stokes_u(aligned_phi, weights=None)
    return aligned_phi, aligned_q, aligned_u


def align_stokes_parameters(q, u, q0, u0):
    """Align the Stokes parameters according to an input polarization model.

    Args
    ----
    q : array_like
        The event-by-event Q Stokes parameters.

    u : array_like
        The event-by-event U Stokes parameters.

    q0 : array_like
        The input model Q Stokes parameters, calculated at the positions of the
        input events.

    u0 : array_like
        The input model U Stokes parameters, calculated at the positions of the
        input events.
    """
    aligned_q = 0.5 * q * q0 + 0.5 * u * u0
    aligned_u = 0.5 * u * q0 - 0.5 * q * u0
    return aligned_q, aligned_u

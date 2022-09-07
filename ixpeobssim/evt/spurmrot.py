#!/usr/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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

"""Spurious modulation correction via the rotation method.
"""

import numpy


# pylint: disable=invalid-name


def delta_phi_ampl(phi, amplitude, phase, harmonic=2):
    """Return the event-by event rotation angle for a given modulation amplitude
    and phase.

    Arguments
    ---------
    phi : array_like
        The original angles.

    amplitude : array_like
        The spurious modulation amplitude.

    phase : array_like
        The spurious modulation phase.

    harmonic : int
        The harmonic to be removed (the default value of 2 corresponds to
        the signal-like spurious modulation).
    """
    return 0.5 * amplitude * (numpy.sin(harmonic * (phi - phase)) + numpy.sin(harmonic * phase))


def delta_phi_stokes(phi, qspur, uspur):
    """Return the event-by event rotation angle for a given spurious Q and U
    Stokes parameters.

    Note this, by definition, does not support harmonics different from 2.

    Arguments
    ---------
    phi : array_like
        The original angles.

    qspur : array_like
        The spurious modulation Q Stokes parameter.

    uspur : array_like
        The spurious modulation U Stokes parameter.
    """
    return 0.5 * (qspur * numpy.sin(2. * phi) + uspur * (1. - numpy.cos(2. * phi)))


def correct_phi_ampl(phi, amplitude, phase, harmonic=2):
    """Azimuthal angle correction in the spurious modulation amplitude and phase
    space.
    """
    return phi + delta_phi_ampl(phi, amplitude, phase, harmonic)


def correct_phi_stokes(phi, qspur, uspur):
    """Azimuthal angle correction in the spurious Stokes parameters space.
    """
    return phi + delta_phi_stokes(phi, qspur, uspur)


def stokes_rotation_angle(q, u, qspur, uspur):
    """Return the rotation angle for the correction in Stokes parameter space.

    Arguments
    ---------
    q : array_like
        The original array of q values.

    u : array_like
        The original array of u values.

    qspur : array_like
        The spurious modulation Q Stokes parameter.

    uspur : array_like
        The spurious modulation U Stokes parameter.
    """
    return 0.5 * u * qspur + (1. - 0.5 * q) * uspur


def correct_stokes_parameters(q, u, qspur, uspur):
    """Correct the spurious modulation via a rotation in Stokes space.
    """
    delta = stokes_rotation_angle(q, u, qspur, uspur)
    c, s = numpy.cos(delta), numpy.sin(delta)
    return q * c - u * s, q * s + u * c

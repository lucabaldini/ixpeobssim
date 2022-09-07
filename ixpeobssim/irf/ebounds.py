#!/usr/bin/env python
#
# Copyright (C) 2018--2021, the ixpeobssim team.
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

from __future__ import print_function, division


import numpy


"""Basic pulse invariant format definition.

All the response functions are defined in true energy between 1 and 12 keV, in
steps of 40 eV.

The pulse invariant spans the 0--15 keV energy band, with the same steps of 40 eV.
(Essentially we are assuming that our response is zero outside the 1--12 keV
energy range, but we do want to capture the effect of the very-low and very-high
energy event that have measured energy below 1 keV or above 12 keV.)
"""

# Definition of the physical space for the response functions---this is the
# ground-truth energy in keV.
ENERGY_MIN = 1.0
ENERGY_MAX = 12.0
ENERGY_STEP = 0.04

# ENERGY_LO and ENERGY_HI arrays for the corresponding columns of the SPECRESP
# extension in the arf files.
ENERGY_LO = numpy.arange(ENERGY_MIN, ENERGY_MAX, ENERGY_STEP)
ENERGY_HI = ENERGY_LO + ENERGY_STEP

# Actual energy binning corresponding to the ENERGY_LO and ENERGY_HI arrays.
ENERGY_BINNING = numpy.append(ENERGY_LO, ENERGY_MAX)

# And, finally, the array containing the bin centers.
ENERGY_GRID = 0.5 * (ENERGY_LO + ENERGY_HI)

# Definition of the channels for the response function.
PI_ENERGY_MIN = 0.0
PI_ENERGY_MAX = 15.0

# Number of channels in the response function
NUM_CHANNELS = int((PI_ENERGY_MAX - PI_ENERGY_MIN) / ENERGY_STEP)
TLMIN = 0
TLMAX = NUM_CHANNELS - 1
PI_BINNING = numpy.linspace(0, NUM_CHANNELS, NUM_CHANNELS + 1)


def energy_to_channel(energy):
    """Convert an energy in keV to the corresponding PI (floating-point) channel.
    """
    return (energy - PI_ENERGY_MIN) / ENERGY_STEP


def digitize_channel(channel, dtype=numpy.int16):
    """Digitize a specific pulse invariant.

    This is essentially rounding an array to the nearest integer and
    cast the result to the specified type.
    """
    return numpy.ndarray.astype(numpy.rint(channel), dtype)


def channel_to_energy(channel):
    """Convert a PI channel to the corresponding physical energy.

    Note we're adding half of the channel step to get the bin center.
    """
    return channel * ENERGY_STEP + PI_ENERGY_MIN + 0.5 * ENERGY_STEP

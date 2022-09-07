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

from __future__ import print_function, division

from ixpeobssim.utils.units_ import mbar_to_atm, celsius_to_kelvin, ZERO_CELSIUS


"""Some relevant physical constants.
"""

# Avogadro number.
N = 6.022140857e23

# Perfect gas constant.
R = 0.08206

# Some atomic masses.
H_MASS = 1.008
HE_MASS = 4.002602
C_MASS = 12.011
O_MASS = 15.999
DME_MASS = 2 * C_MASS + 6 * H_MASS + O_MASS

# Some densities (in cm^{-3}).
BE_DENSITY = 1.85
AL_DENSITY = 2.70

# Density of gaseous DME at 1 atm, 0 deg C in g/cm3
# See https://bitbucket.org/ixpesw/ixpeobssim/issues/284/
DME_REF_DENSITY = 2.1146e-3


def perf_gas_molar_volume(temperature, pressure):
    """Return the molar volume (in cm^3) at a given temperature and pressure.

    This is using the law of perfect gases.

    Arguments
    ---------
    temperature : float
        The temperature in degrees C.

    pressure : float
        The pressure in mbar.
    """
    return R * celsius_to_kelvin(temperature) / mbar_to_atm(pressure) * 1000.


def perf_gas_density(mass, temperature, pressure):
    """Return the density (in g cm^{-3}) for a perfect gas of a given
    molecular mass at a given temperature and pressure.

    Arguments
    ---------
    mass : float
        The molecular mass of the gas.

    temperature : float
        The temperature in degrees C.

    pressure : float
        The pressure in mbar.
    """
    return mass / perf_gas_molar_volume(temperature, pressure)


def _dme_density_scaling(temperature, pressure, index):
    """Generic function for the DME density scaling.

    This assumes a power-law parametrization, giving the standard law of
    perfect gases when index = 1.

    Arguments
    ---------
    temperature : float
        The temperature in degrees C.

    pressure : float
        The pressure in mbar.
    """
    return DME_REF_DENSITY * mbar_to_atm(pressure) * \
        (ZERO_CELSIUS / celsius_to_kelvin(temperature))**index


def dme_density_perfect(temperature, pressure):
    """Return the density for the gaseous DME as a function of the temperature
    and pressure starting following the law of perfect gases.

    Arguments
    ---------
    temperature : float
        The temperature in degrees C.

    pressure : float
        The pressure in mbar.
    """
    return _dme_density_scaling(temperature, pressure, 1.)


def dme_density(temperature, pressure):
    """Return the density for the gaseous DME as a function of the temperature
    and pressure based on our custom parametrization.

    See https://bitbucket.org/ixpesw/ixpeobssim/issues/284/ for more details.

    Arguments
    ---------
    temperature : float
        The temperature in degrees C.

    pressure : float
        The pressure in mbar.
    """
    return _dme_density_scaling(temperature, pressure, 1.10042278)

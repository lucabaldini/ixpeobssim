#!/urs/bin/env python
#
# Copyright (C) 2016, the ixpeobssim team.
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

"""Facilities for conversions across different units.
"""

from __future__ import print_function, division


FACTOR_KEV_TO_ERG = 6.2415096471204e8
FACTOR_ERGCMS_TO_MCRAB = 0.5e11
FACTOR_MBAR_TO_ATM = 1013.249977
ZERO_CELSIUS = 273.15


def erg_to_keV(val):
    """Convert erg to keV.
    """
    return val * FACTOR_KEV_TO_ERG

def keV_to_erg(val):
    """Convert keV to erg.
    """
    return val / FACTOR_KEV_TO_ERG

def ergcms_to_mcrab(val):
    """Convert a flux in erg per cm square into mcrab.
    """
    return val * FACTOR_ERGCMS_TO_MCRAB

def arcmin_to_degrees(val):
    """Convert arcminutes to degrees.
    """
    return val / 60.

def degrees_to_arcmin(val):
    """Convert degrees to arcminutes.
    """
    return val * 60.

def arcsec_to_degrees(val):
    """Convert arcseconds to degrees.
    """
    return val / 3600.

def degrees_to_arcsec(val):
    """Convert degrees to arcseconds.
    """
    return val * 3600.

def arcsec_to_arcmin(val):
    """Convert arcseconds to arcminutes.
    """
    return val / 60.

def arcmin_to_arcsec(val):
    """Convert arcminutes to arcseconds.
    """
    return val * 60.

def mbar_to_atm(val):
    """Convert mbar to atm.
    """
    return val / FACTOR_MBAR_TO_ATM

def atm_to_mbar(val):
    """Convert atm to mbar.
    """
    return val * FACTOR_MBAR_TO_ATM

def celsius_to_kelvin(val):
    """Convert Celsius degrees to Kelvin degrees.
    """
    return val + ZERO_CELSIUS

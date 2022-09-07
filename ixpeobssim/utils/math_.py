#!/usr/bin/env python
#
# Copyright (C) 2015, the ixpeobssim team.
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
import sys


def fold_angle_deg(phi):
    """Fold an azimuthal angle (in degrees) within [-180, 180] to [-90, 90].

    Parameters
    ----------
    phi : array_like or scalar
        The input angle or array of values.
    """
    return phi +\
        180. * numpy.logical_and(phi >= -180, phi < -90.) -\
        180. * numpy.logical_and(phi > 90., phi <= 180.)


def fold_angle_rad(phi):
    """Fold an azimuthal angle (in degrees) within [-pi, pi] to [-pi/2, pi/2].

    Parameters
    ----------
    phi : array_like or scalar
        The input angle or array of values.
    """
    return phi +\
        numpy.pi * numpy.logical_and(phi >= -numpy.pi, phi < -0.5 * numpy.pi) -\
        numpy.pi * numpy.logical_and(phi > 0.5 * numpy.pi, phi <= numpy.pi)


def modulo_2pi(phi):
    """Compute the modulo operation bringing the output angles (in radians)
    into the interval [-pi, pi].

    Parameters
    ----------
    phi : array_like or scalar
        The input angle or array of values.
    """
    return phi + 2 * numpy.pi * (phi < -numpy.pi) - 2 * numpy.pi * (phi > numpy.pi)


def format_value(value, precision=3):
    """Format a number with a reasonable precision
    """
    if isinstance(value, str):
        return value
    else:
        fmt = '%%.%dg' % precision
        return fmt % value


def decimal_places(val):
    """Calculate the number of decimal places so that a given value is rounded
    to exactly two signficant digits.

    Note that we add epsilon to the argument of the logarithm in such a way
    that, e.g., 0.001 is converted to 0.0010 and not 0.00100. For values greater
    than 99 this number is negative.
    """
    return 1 - int(numpy.log10(val + sys.float_info.epsilon)) + 1 * (val < 1.)


def decimal_power(val):
    """Calculate the order of magnitude of a given value,i.e., the largest
    power of ten smaller than the value.
    """
    return int(numpy.log10(val + sys.float_info.epsilon)) - 1 * (val < 1.)


def format_value_error(value, error, pm='+/-', max_dec_places=6):
    """Format a measurement with the proper number of significant digits.
    """
    value = float(value)
    error = float(error)
    assert error >= 0
    if error == 0:
        return '%s %s 0' % (format_value(value), pm)
    if error == numpy.inf:
        return '%s %s inf' % (format_value(value), pm)
    dec_places = decimal_places(error)
    if dec_places >= 0 and dec_places <= max_dec_places:
        fmt = '%%.%df %s %%.%df' % (dec_places, pm, dec_places)
    else:
        p = decimal_power(abs(value))
        scale = 10 ** p
        value /= scale
        error /= scale
        dec_places = decimal_places(error)
        if dec_places > 0:
            if p > 0:
                exp = 'e+%02d' % p
            else:
                exp = 'e-%02d' % abs(p)
            fmt = '%%.%df%s %s %%.%df%s' %\
                  (dec_places, exp, pm, dec_places, exp)
        else:
            fmt = '%%d %s %%d' % pm
    return fmt % (value, error)


def weighted_average(values, weights=None):
    """Return the weighted average of an array of values.

    This is simply a wrapper over the numpy.average() function.
    """
    return numpy.average(values, weights=weights)

#!/urs/bin/env python
#
# Copyright (C) 2019, the ixpeobssim team.
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

"""Spacecraft-related facilities.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.core.spline import xInterpolatedUnivariateSpline


def period_to_omega(period):
    """Convert a characteristic period (in s) into the corresponding omega
    (in rad/s).
    """
    assert period > 0
    return 2. * numpy.pi / period

def dithering_pattern(amplitude=1.6, period_a=907., period_x=101., period_y=449.):
    """Implementation of the full dithering pattern as per the report by
    Allyn Tennant and Kurt Dietz linked in the issue page:
    https://bitbucket.org/ixpesw/ixpeobssim/issues/193

    Note this returns an anonymous function that can be evaluated into a generic
    array of time values, the basic usage being:

    >>> t = numpy.linspace(0., 10000., 10001)
    >>> dithering = dithering_pattern()
    >>> x, y = dithering(t)

    Arguments
    ---------
    amplitude : float
        The dithering amplitude in arcminutes.

    period_a : float
        The main dither period.

    period_x : float
        The x dithering period.

    period_y : float
        The y dithering period.
    """
    omega_a = period_to_omega(period_a)
    omega_x = period_to_omega(period_x)
    omega_y = period_to_omega(period_y)
    x = lambda t: amplitude * numpy.cos(omega_a * t) * numpy.cos(omega_x * t)
    y = lambda t: amplitude * numpy.sin(omega_a * t) * numpy.sin(omega_y * t)
    return lambda t: (x(t), y(t))


def triangular_wave(x, amplitude, period):
    """Basic description of a (symmetric) triangular wave with a given period
    and amplitude.

    The basic expression is taken from
    https://en.wikipedia.org/wiki/Triangular_distribution
    and, in this form, the function evaluates to 0 for x = 0.
    """
    half_period = 0.5 * period
    # Take the input x array modulo the period.
    x = numpy.mod(x, period)
    # Use the wikipedia formula.
    return amplitude * (half_period - numpy.abs(half_period - x)) / half_period


def pow_triangular_wave(x, amplitude, period, exponent=0.5):
    """Modified triangual wave, where the relative values are raised to a
    generic exponent, in such a way that the values of the maxima and minima are
    preserved.
    """
    return amplitude * (triangular_wave(x, amplitude, period) / amplitude) ** exponent


def spiral_dithering_pattern(amplitude=1.6, period_theta=100., period_r=970.,
                             exponent=0.5):
    """Alternative, spiral-like dithering pattern.
    """
    omega = period_to_omega(period_theta)
    r = lambda t: pow_triangular_wave(t, amplitude, period_r, exponent)
    x = lambda t: r(t) * numpy.cos(omega * t)
    y = lambda t: r(t) * numpy.sin(omega * t)
    return lambda t: (x(t), y(t))


def pointing_splines(sc_data):
    """Build a pair of R.A. and Dec dithering splines starting from as set
    of spacecraft data.

    This can be used to recover the actual pointing as a function of time, given
    a SC_DATA binary table.
    """
    met = sc_data['MET']
    ra = sc_data['RA_PNT']
    dec = sc_data['DEC_PNT']
    ra_spline = xInterpolatedUnivariateSpline(met, ra, xlabel='MET [s]', ylabel='Ditehring R.A.')
    dec_spline = xInterpolatedUnivariateSpline(met, dec, xlabel='MET [s]', ylabel='Ditehring Dec')
    return ra_spline, dec_spline


def pointing_direction(sc_data, met):
    """Return the pointing direction at a given array of MET.
    """
    ra_spline, dec_spline = pointing_splines(sc_data)
    return ra_spline(met), dec_spline(met)

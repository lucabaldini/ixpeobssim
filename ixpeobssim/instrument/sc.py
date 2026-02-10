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

from astropy.io import fits

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


def mask_to_gti(time, mask, *, t_start=None, t_end=None):
    """
    Convert a boolean mask into half-open GTIs [t_start, t_stop).

    If t_start or t_end are provided, they define the observation
    boundaries. Otherwise, time[0] and time[-1] are used.

    Arguments
    ----------
    time : ndarray
        1D monotonically increasing time array.
    mask : ndarray of bool
        True inside GTIs.
    t_start : float or None, optional
        Observation start time. Must satisfy t_start <= time[0].
    t_end : float or None, optional
        Observation end time. Must satisfy t_end >= time[-1].

    Returns
    -------
    tstart, tstop : ndarray
        Start and stop times of GTIs (half-open).
    """
    assert mask.ndim == 1
    assert time.size == mask.size
    if t_start is not None and t_start > time[0]:
        raise ValueError(
            't_start must be <= time[0] '
            f'(t_start={t_start}, time[0]={time[0]})'
        )
    if t_end is not None and t_end < time[-1]:
        raise ValueError(
            't_end must be >= time[-1] '
            f'(t_end={t_end}, time[-1]={time[-1]})'
        )
    # Convert to int and perform edge detection
    mask_i = mask.astype(int)
    edges = numpy.ediff1d(mask_i)
    # Indices where a GTI starts (0 → 1 transition)
    start_idx = numpy.where(edges > 0)[0] + 1
    # Indices where a GTI ends (1 → 0 transition)
    stop_idx  = numpy.where(edges < 0)[0] + 1

    # Handle open interval at start
    if mask_i[0]:
        start_idx = numpy.insert(start_idx, 0, -1)
    # Handle open interval at end
    if mask_i[-1]:
        stop_idx = numpy.append(stop_idx, mask_i.size)
    assert start_idx.size == stop_idx.size
    # Map indices to times
    tstart = numpy.empty_like(start_idx, dtype=time.dtype)
    tstop  = numpy.empty_like(stop_idx,  dtype=time.dtype)

    m_start = (start_idx == -1)
    m_stop  = (stop_idx == mask_i.size)

    tstart[m_start] = time[0] if t_start is None else t_start
    tstart[~m_start] = time[start_idx[~m_start]]

    tstop[m_stop] = time[-1] if t_end is None else t_end
    tstop[~m_stop] = time[stop_idx[~m_stop]]

    return tstart, tstop


def create_ineclipse_gtis(*hk_file_paths, complement=False):
    """
    Create arrays of start and stop times defining Good Time Intervals (GTIs)
    corresponding to spacecraft eclipse conditions.
        
    By default (complement=False), the function returns intervals when the
    spacecraft is not affected by the sun.

    If complement=True, the returned intervals correspond to everything else.
    """
    # Lists collecting GTI start/stop times from all input HK files
    tstart_all = []
    tstop_all = []
    for file_path in hk_file_paths:
        with fits.open(file_path) as hdul:
            data = hdul['HK'].data
            # Angular distances (sign encodes geometry / visibility)
            adsec2ecl = data['ADSEC2ECL']
            adsec2sun = data['ADSEC2SUN']
            time = data['TIME']
            if not numpy.all(numpy.diff(time) >= 0):
                raise ValueError(f'TIME column in {file_path} is not monotonic')
            # Build a boolean mask selecting the desired condition:
            #  - complement=False → eclipse / not INSUN
            #  - complement=True  → INSUN
            if complement:
                mask = (adsec2ecl < 0) | (adsec2sun >= 0)
            else:
                mask = (adsec2sun < 0) * (adsec2ecl >= 0)
            tstart, tstop = mask_to_gti(time, mask)
            tstart_all.append(tstart)
            tstop_all.append(tstop)
    # Concatenate GTIs from all input files into single arrays
    return numpy.concatenate(tstart_all), numpy.concatenate(tstop_all)

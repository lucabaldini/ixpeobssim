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

"""
Module containing functions and classes to manage time delays
"""

from __future__ import print_function, division

import numpy
import astropy.constants as const

from ixpeobssim.srcmodel.ephemeris import get_earth_barycentric_ephemeris, get_eccentric_anomaly
from ixpeobssim.utils.time_ import met_to_mjd, mjd_to_met
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline

# pylint: disable=invalid-name, too-many-arguments, too-many-locals

"""
Time delays in the Solar System
"""

def roemers_f(t_, **kwargs):
    """Roemer delays formula referred to Solar System Barycenter.

    Warning
    -------
    IXPE satellite ephemeris not implemented yet.

    Args
    ----
    t_ : float
       Terrestrial time in MET

    **kwargs :
       ra : Right ascension of the source
       dec : Declination of the source

    Returns
    -------
    array float
        Return the earth barycentric vector
    """

    earth_vector = get_earth_barycentric_ephemeris(t_)

    c_ = const.c.value # vacuum speed of light
    # source coordinates
    lat = numpy.deg2rad(kwargs['ra'])
    lon = numpy.deg2rad(kwargs['dec'])

    # source point unit vector from SSB
    s_ = numpy.array([numpy.cos(lat)*numpy.cos(lon),
                      numpy.cos(lat)*numpy.sin(lon),
                      numpy.sin(lat)])

    # delay is the scalar product between earth SSB vector and source position
    func = - (1. / c_) * numpy.dot(earth_vector, s_)
    return func

def einsteins_f(t_, **kwargs):
    """Einstein delays formula referred to Solar System Barycenter.

    Warning
    -------
    First implementation of satellite location, temporarily avoided.
    """

    # Here we put an order of magnitude for IXPE orbit radius
    #tt_ = Time(met_to_mjd(t_), format='mjd', scale='tt',
    #           location=(0.*u.deg, 0.*u.deg, 5.4e5*u.meter))
    #tdb_ = mjd_to_met(tt_.tdb.value)
    #func = tdb_ - t_
    func = 0.
    return func

def shapiros_f(t_, **kwargs):
    """Shapiro delays formula referred to Solar System Barycenter.

    Warning
    -------
    Not implemented yet.
    """
    func = 0.
    return func

def alls_f(t_, **kwargs):
    """All delays referred to Solar System Barycenter.
    """
    func = roemers_f(t_, **kwargs) + \
           einsteins_f(t_, **kwargs) + \
           shapiros_f(t_, **kwargs)
    return func


"""
Time delays in the Binary System
"""

def roemerb_f(t_, eph, **kwargs):
    """Roemer delays formula referred to Binary System.

    Warning
    -------
    IXPE satellite ephemeris not implemented yet.

    Args
    ----
    t_ : float
       Terrestrial time in MET

    eph : xOrbitalEphemeris object

    **kwargs :
       ra : Right ascension of the source
       dec : Declination of the source

    Returns
    -------
    array float
        Return the earth barycentric vector
    """

    if eph.model == 'ELL1':
        # Roemer delay argument phi
        dt = (t_ - eph.t_orbital)
        phi = 2. * numpy.pi * eph.omega_orb(dt) * dt
        # Here, I use the ELL1 model to apply Roemer delay to the TOAs/times vector
        func = eph.axis_proj * (numpy.sin(phi) + \
               eph.eps1 * 0.5 * numpy.sin(2 * phi) - \
               eph.eps2 * 0.5 * numpy.cos(2 * phi))
        return func
    # Here we solve Kepler Equation for Roemer binary delays
    # Eccentric anomaly from K.E. Solution
    U_ = get_eccentric_anomaly(t_, eph)
    # Longitude of periastron in radians
    OM_ = numpy.deg2rad(eph.lon_periast)
    func = eph.axis_proj * (numpy.cos(U_) - eph.ecc) * numpy.sin(OM_) + \
           eph.axis_proj * numpy.sin(U_) * numpy.sqrt(1 - eph.ecc**2) * numpy.sin(OM_)
    return func

def einsteinb_f(t_, **kwargs):
    """Einstein delays formula referred to Binary System.

    Warning
    -------
    Not implemented yet.
    """
    func = 0.
    return func

def shapirob_f(t_, **kwargs):
    """Shapiro delays formula referred to Binary System.

    Warning
    -------
    Not implemented yet.
    """
    func = 0.
    return func

def allb_f(t_, eph, **kwargs):
    """All delays referred to Binary System.
    """
    func = roemerb_f(t_, eph, **kwargs) + \
           einsteinb_f(t_, **kwargs) + \
           shapirob_f(t_, **kwargs)
    return func

"""
Dictionary for delays in Solar System (SSB)
"""

SSB = {'roemer': roemers_f,
       'einstein': einsteins_f,
       'shapiro': shapiros_f,
       'all': alls_f}
"""
Dictionary for delays in Binary System (BS)
"""

BS = {'roemer': roemerb_f,
      'einstein': einsteinb_f,
      'shapiro': shapirob_f,
      'all': allb_f}



class xTDelays:
    """Convenience class encapsulating times handling.

    Default time in input is in mjd unit, both mjd and met time scale are managed.

    Args
    ----
    mjdvalue : array float
       array of time values in MJD unit

    metvalue : array float
       array of time values in MET unit

    unit : string
       time unit in MET or MJD

    name : string
       name of the model
    """

    def __init__(self, mjdvalue, name='', unit='MJD'):
        """Constructor.
        """

        assert unit.upper() in ('MET', 'MJD')

        self.mjdvalue = numpy.sort(mjdvalue)
        self.metvalue = mjd_to_met(self.mjdvalue)
        if unit in ('met', 'MET'):
            self.metvalue = self.mjdvalue
            self.mjdvalue = met_to_mjd(self.mjdvalue)
        self.unit = unit
        self.name = name

    def mjd_min(self):
        """Return the minimum in mjd unit.
        """
        return self.mjdvalue[0]

    def mjd_max(self):
        """Return the maximum in mjd unit.
        """
        return self.mjdvalue[-1]

    def met_min(self):
        """Return the minimum in met unit.
        """
        return self.metvalue[0]

    def met_max(self):
        """Return the maximum in met unit.
        """
        return self.metvalue[-1]

    def __len__(self):
        """Return the lenght of the array.
        """
        return len(self.mjdvalue)

    def __str__(self):
        """String formatting.
        """
        return 'xTDelays = %s s, ref = %s, name = %s' % (self.mjdvalue, self.unit, self.name)

    def apply_delay(self, ephemeris, ra, dec, delay='all', binary=False):
        """Return barycentered times, given an ephemeris and the source coordinates.

        Args
        ----
        ephemeris : xOrbitalEphemeris object

        ra : float
           Right ascension of the source

        dec : float
           Declination of the source

        delay : string
           define the delay function to use: roemer, einstein, shapiro or all

        binary : Boolean
           When true time delay is applied in the binary system, when false in the solar system
        """
        delta = SSB[delay](self.metvalue, ra=ra, dec=dec)
        if binary:
            delta += BS[delay](self.metvalue, ephemeris, ra=ra, dec=dec)
        self.__init__(self.metvalue + delta, unit='met', name='%s_DELAY' % self.name)
        return self.metvalue

    def apply_decorr(self, ephemeris, ra, dec, samples=500, delay='all', binary=False):
        """Return delay effected times.

        Args
        ----
        samples : int
           Number of samples per period
        """
        # We need enough samples for each periodicity
        if hasattr(ephemeris, 'model'):
            num_periods = (self.met_max() - self.met_min()) / ephemeris.p_orbital
            if num_periods > 1:
                samples = int(samples * num_periods)
        t_ = numpy.linspace(self.met_min(), self.met_max() * 1.1, samples)
        delta = xInterpolatedUnivariateSpline(t_, t_)
        delta += xInterpolatedUnivariateSpline(t_, SSB[delay](t_, ra=ra, dec=dec))
        if binary:
            delta += xInterpolatedUnivariateSpline(t_, BS[delay](t_, ephemeris, ra=ra, dec=dec))
        metvalue = delta.inverse()(self.metvalue)
        self.__init__(metvalue, unit='met', name='%s_DECORR' % self.name)
        return self.metvalue

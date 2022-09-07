# Copyright (C) 2015--2022, the ixpeobssim team.
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

"""Module containing ephemeris utilities, classes and functions
"""

from __future__ import print_function, division

import numpy
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric
from astropy.time import Time

from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.utils.time_ import met_to_mjd, mjd_to_met, days_to_seconds


# pylint: disable=invalid-name, too-many-arguments, too-many-locals, no-member
# pylint: disable=too-many-instance-attributes


def get_earth_barycentric_ephemeris(met):
    """Return earth barycentric vector given a set of terrestrial times
    referred to MET in seconds.

    Args
    ----
    met : array_like
       Array of terrestrial times in MET

    Returns
    -------
    array
        Return the earth barycentric vector
    """
    #Set ephemeris
    solar_system_ephemeris.set('builtin') #ERFA library
    Tt = Time(met_to_mjd(met), format='mjd', scale='tt')
    #Barycentric (ICRS) position of the body in cartesian coordinates
    earth_b = get_body_barycentric('earth', Tt)
    #earth vector
    e_v = earth_b.get_xyz()
    e_v = e_v.to(u.meter).value
    earth_vector = numpy.matrix.transpose(numpy.array(e_v))
    return earth_vector


def get_eccentric_anomaly(t_, orb_eph, periods=100, samples=10, plot=False):
    """At different time of emission of pulsar we calculate eccentric anomaly of binary system.
       Eccentric Anomaly in radiants.

    Args
    ----
    t_ : array float
       Array of times at which evaluate the eccentric anomaly

    orb_eph : class xOrbitalEphemeris
       Orbital ephemeris

    periods : int
       Number of periods to use to evaluate the eccentric anomaly

    samples : int
       Number of bins in the eccentric anomaly evaluation

    plot : Boolean
       If True plot the inverse function

    Returns
    -------
    class spline
       returns the eccentric anomaly as a function of pulsar emission time
    """
    # EA sampling
    EA_ = numpy.linspace(0, periods * 2 * numpy.pi, samples)
    # Angular velocity, assumed almost the same for all over observation
    n_ = 2 * numpy.pi * orb_eph.omega_orb(t_.min())
    # We solve the inverse of toes function, using spline method
    toes_ = (EA_ - orb_eph.ecc * numpy.sin(EA_)) / n_ + orb_eph.t_orbital
    fmt = dict(xlabel='Ecc Anomaly [rad]', ylabel='TOES [s]')
    toes_splinefunc = xInterpolatedUnivariateSpline(EA_, toes_, **fmt)
    EA_splinefunc = toes_splinefunc.inverse()
    if plot:
        toes_splinefunc.plot()
        plt.figure()
        EA_splinefunc.plot()
    return EA_splinefunc(t_)


def read_par_file(parfile_path):
    """Method to read a parameter file and to get pulsar parameters.

    The following parameters are accepted:

    * PEPOCH   Epoch of period/frequency (MJD)
    * POSEPOCH Epoch of position measuremet (MJD)
    * F0       Pulsar rotation frequency (s^-1)
    * F1       Pulsar rotation frequency derivative (s^-2)
    * F2       Pulsar rotation frequency second derivative (s^-3)
    * P0       Pulsar period (s).
    * P1       Pulsar period derivative (10^-15).
    * DM       Dispersion measure (pc cm^-3)
    * A1       Projected pulsar semi-major axis of orbit (lt-s)
    * E/ECC    Eccentricity of orbit
    * T0       Epoch of periastron passage of orbit (MJD)
    * TASC     Epoch of ascending node passage (MJD)
    * PB       Period of orbit (days)
    * OM       Longitude of periastron passage (deg)
    * EPS1     First Laplace parameter [E x sin(OM)]
    * EPS2     Second Laplace parameter [E x cos(OM)]
    * EPS1DOT  Time derivative of EPS1
    * EPS2DOT  Time derivative of EPS2
    * OMDOT    Rate of periastron advance (deg/yr)
    * PBDOT    Rate of change of orbital period
    * XDOT     Rate of change of projected semi-major axis
    * ECCDOT   Rate of change of eccentricity

    For further information refer to the TEMPO2 user manual.

    .. note::

       Parameter files are commonly used in pulsar timing analysis, because they
       contain precise information about pulsar position, spin-down rate, orbital
       parameters (if in binary systems), etc. The structure matches with a
       dictionary fill-in method.

    Args
    ----
    parfile_path : str
        Path to the .par file containing the ephemeris

    Returns
    -------
    dictionary
        Return a dictionary containing the pulsar parameters
    """
    assert parfile_path.endswith('.par')
    parameters = {}
    with open(parfile_path) as input_file:
        rows = input_file.readlines()
    for line in [row.strip() for row in rows]:
        # Skip comments
        if line[0] == '#':
            continue
        elements = line.split()
        # Skip blank lines
        if len(elements) == 0:
            continue
        _key = elements[0]
        if len(elements) == 2:
            try:
                parameters[_key] = float(elements[1])
            except ValueError:
                parameters[_key] = elements[1]
        # Some lines have errors
        if len(elements) == 3:
            if elements[2] not in ['0', '1']:
                try:
                    parameters[_key] = float(elements[1])
                    parameters[_key + '_ERR'] = float(elements[2])
                except ValueError:
                    pass
            else:
                try:
                    parameters[_key] = float(elements[1])
                except ValueError:
                    pass
        if len(elements) == 4:
            try:
                parameters[_key] = float(elements[1])
                parameters[_key + '_ERR'] = float(elements[3])
            except ValueError:
                pass
    return parameters



def time_list(pulse_shape, start_met, ephemeris, n_events, duration):
    """Time events extracted with the acceptance rejection method.

    Args
    ----

    pulse_shape :
       Pulse profile shape from pulse_pol_from_harmonics_spline in ixpeobssim.srcmodel.polarization

    start_met : float
       Starting time

    ephemeris : xOrbitalEphemeris object
        The orbital ehemeris to be used.

    n_events : integer
       Number of events to simulate

    duration : float
       Time of observation
    """
    times = numpy.random.uniform(start_met, start_met + duration, n_events * 3)
    times = numpy.sort(times)
    # Timing model applied to the times
    ph = phase_function(times, start_met, ephemeris.nu(start_met), ephemeris.nudot(start_met))
    ph = ph - numpy.floor(ph)
    max_ = numpy.amax(pulse_shape(ph))

    # Randomly flat distributed amplitude for pulse shape
    amplitudes = numpy.random.uniform(0., max_, n_events * 3)

    # Comparison check and selection of good amplitude values
    good = amplitudes < pulse_shape(ph)
    return times[good][:n_events]



class xEphemeris:

    """Convenience class encapsulating the ephemeris for a periodic source.

    This is meant to capture effects up to the second time derivative of the
    frequency.

    The basic task of the class is to handle the conversion from mission elapsed
    time to phase and vice versa. This is performed through the met_to_phase(),
    and phase_to_met() interfaces. Note that, in this context, by `phase` we
    mean a free-running quantity, while by `pulse phase` we mean the fractional
    part of the former, in the [0, 1[ interval.

    The direct transformation (met to phase or pulse phase) is performed through
    a simple taylor expansion

    .. math::

       \\phi = \\nu_0 (t - t_0) + \\frac{1}{2} \\dot\\nu_0 (t - t_0)^2 +
       \\frac{1}{6} \\ddot\\nu (t - t_0)^3

    while the inverse (phase to met) is done through an interpolating spline.

    In addition, given a pulse shape in the form of a xUnivariateGenerator, the
    class is able to extract a set of random arrays of met and pulse phase
    via the rvs() method. This is used in ixpeobssim.srcmodel.roi for
    simulating periodic sources.

    Args
    ----
    met0 : float
        The reference mission elapsed time.

    nu0 : float
        Rotational frequency at met0.

    nudot0 : float
        First detivative of the rotational frequency at met0.

    nuddot : float
        Second derivative of the rotational frequency (constant).
    """

    def __init__(self, met0, nu0, nudot0=0., nuddot=0.):
        """Constructor.
        """
        self.met0 = met0
        self.nu0 = nu0
        self.nudot0 = nudot0
        self.nuddot = nuddot

    def dict(self):
        """Return the ephemeris content in the form of a dictionary.

        This is useful, e.g., when calling xpphase passing an ephemeris object
        to set the fields programmatically. The reason why we're wrapping this
        and not using __dict__ directly is to be able to accommodate possible
        future changes in either the class layout or the options for the
        external programs.
        """
        return self.__dict__

    def _dt(self, met):
        """Return the time difference wrt the reference met.
        """
        return met - self.met0

    def nu(self, met):
        """Return the source frequency at a given time.
        """
        dt = self._dt(met)
        return self.nu0 + self.nudot0 * dt + self.nuddot * (dt**2.) / 2.

    def nudot(self, met):
        """Return the source first derivative of frequency at a given time.
        """
        dt = self._dt(met)
        return self.nudot0 + self.nuddot * dt

    def period(self, met):
        """Return the source period at a given time.
        """
        return 1. / self.nu(met)

    def met_to_phase(self, met):
        """Return the phase for a given value of array of MET values.

        Note the phase, here, is free-running, i.e. not bound between 0 and 1.
        For a simple ephemeris with zero time derivatives this is returning
        a straight line with a slope equal to the (constant) frequency.
        """
        dt = self._dt(met)
        return self.nu0 * dt + self.nudot0 * (dt**2.) / 2. + self.nuddot * (dt**3.) / 6.

    def fold(self, met, start_met, phi0=0.):
        """Fold the MET values to return the corresponding arrays of pulse phase.

        .. warning::

            This is creating a new xEphemeris object with the reference MET set
            to start_met, and the frequency derivatives changed accordingly, and
            then calling the met_to_phase() method of this new object. This is
            supposed to roundtrip with the rvs() method to a decent accuracy,
            and is therefore the interface used in the xpphase application.
            Note, however, that there is a lot of numerical (more or less random)
            rounding going on under the hood, and this is not intended
            for "professional" use. (In a sense, the function internals are
            tweaked specifically for the numerical noise to play in such a way
            that we're making approximately the same errors in rvs() and fold()).
        """
        met0 = start_met
        nu0 = self.nu(start_met)
        nudot0 = self.nudot(start_met)
        phase = xEphemeris(met0, nu0, nudot0, self.nuddot).met_to_phase(met)
        pulse_phase = numpy.mod(phase + phi0, 1)
        return pulse_phase

    def phase_spline(self, start_met, duration, num_points=1000):
        """Return an interpolated univariate spline with the phase values for
        a given time range.

        This is mainly used in the inverted fashion to convert from phase to met.
        """
        met = numpy.linspace(start_met, start_met + duration, num_points)
        fmt = dict(xlabel='Time [s]', ylabel='Phase')
        return xInterpolatedUnivariateSpline(met, self.met_to_phase(met), **fmt)

    def phase_spline_inverse(self, start_met, duration, num_points=1000):
        """Return the inverse of the phase spline for a given time range.

        This can effectively be used for calculating the proper mission elapsed
        times starting from a bunch of phase values.
        """
        return self.phase_spline(start_met, duration, num_points).inverse()

    def _zero_order_met_guess(self, phase):
        """Provide a zero-order guess of the mission elapsed time(s) corresponding
        to a particolar value  or array of phases.
        """
        return phase / self.nu0 + self.met0

    def phase_to_met(self, phase, start_met=None, duration=None, num_points=1000):
        """Return the mission elapsed times for given value or array of phases.

        This is essentially the inverse of phase(), and the two are supposed to
        roundtrip with each other. The reason why this function is slightly more
        complicated is that the former cannot be inverted analytically, and
        since the source frequency changes with time, it is not trivial to map
        directly the minimum and maximum phase to the corresponding met values,
        which is needed to create the interpolating spline for the conversion.

        The start_met and stop_met arguments are provided to control the spline
        operating the conversion, and should be provided if available---in this
        case the span of the spline is optimal and it is guaranteed that no
        extrapolation is needed. Both arguments default to None, in which case a
        zero-order guess is made as to what their value shoule be.

        Args
        ----
        phase : array_like
            The phase value(s) to be converted in met.

        start_met : float
            The starting met value for the spline operating the conversion.

        duration : float
            The overall span for the spline operating the conversion.

        num_points : int
            The number of points for spline operating the conversion.
        """
        if start_met is None:
            start_met = self._zero_order_met_guess(phase.min())
        if duration is None:
            duration = self._zero_order_met_guess(phase.max()) - start_met
        return self.phase_spline_inverse(start_met, duration, num_points)(phase)

    def rvs(self, pulse_profile, start_met, duration, num_events):
        """Return two (sorted) arrays of pulse phase and met for a given
        pulse profile, observation interval and number of events.

        Args
        ----
        pulse_profile : xUnivariateGenerator instance
            The pulse target profile.

        start_met : float
            The initial met for the observation.

        duration : float
            The duration of the observation.

        num_events : float
            The total number of events to be simulated.
        """
        # Calculate the number of periods spanned by the target met interval.
        # Note that we split the calculation into the integer and fractional
        # parts, that need to be handled separately.
        phase_min = self.met_to_phase(start_met)
        phase_max = self.met_to_phase(start_met + duration)
        delta_phase = phase_max - phase_min
        reminder, num_periods = numpy.modf(delta_phase)
        # Calculate the number of events to be extracted in the last fraction
        # of a period.
        num_extra_events = numpy.round(num_events * reminder / (delta_phase)).astype(int)
        # Step 1: handle the first n full periods: extract two random numbers,
        # one for the (integral) number of periods and one for the (fractional)
        # pulse phase.
        period = numpy.random.randint(0, num_periods, num_events - num_extra_events)
        pulse_phase = pulse_profile.rvs(period.shape)
        phase = phase_min + period + pulse_phase
        # Step 2: if the target duration does not span exactly an integer number
        # of periods, we take care of the last fractional part.
        if num_extra_events > 0:
            _pulse_phase = pulse_profile.rvs_bounded(num_extra_events, rvmax=reminder)
            pulse_phase = numpy.append(pulse_phase, _pulse_phase)
            phase = numpy.append(phase, phase_min + num_periods + _pulse_phase)
        # Transfor the phase to MET values.
        met = self.phase_to_met(phase, start_met, duration)
        # Sort the event times (and phases accordingly).
        mask = numpy.argsort(met)
        pulse_phase = pulse_phase[mask]
        met = met[mask]
        # And we're good to go!
        return pulse_phase, met

    @classmethod
    def from_file(cls, parfile_path):
        """Read a set of attributes given a parameter file.
        """
        par_dict = read_par_file(parfile_path)
        params = {}
        params['met0'] = mjd_to_met(par_dict['PEPOCH'])
        params['nu0'] = par_dict['F0']
        params['nudot0'] = par_dict.get('F1', 0.)
        params['nuddot'] = par_dict.get('F2', 0.)
        return cls(**params)

    def __str__(self):
        """String formatting.
        """
        return '%s %s' % (self.__class__.__name__, self.__dict__)



def phase_function(met, met0, nu0, nudot0=0., nuddot=0.):
    """Return the rotational phase for a given observational time.

    Warning
    -------
    This function is deprecated in favor of the xEphemeris.met_to_phase() class
    method, and is (temporarily) kept for backward compatibility.
    """
    return xEphemeris(met0, nu0, nudot0, nuddot).met_to_phase(met)



class xOrbitalEphemeris(xEphemeris):

    """Convenience class encapsulating an orbital ephemeris.

    Arguments
    ---------
    met0 : float
        Epoch of period/frequency (MET)

    nu0 : float
        Pulsar rotation frequency (s^-1)

    t_orbital : float
        Epoch of ascending node passage or Epoch of periastron passage of orbit (MET)

    p_orbital : float
        Period of orbit (s)

    axis_proj : float
        Projected pulsar semi-major axis of orbit (lt-s)

    ecc : float
        Eccentricity of orbit

    lon_periast : float
        Longitude of periastron passage (deg)

    nudot0 : float
        Pulsar rotation frequency derivative (s^-2)

    nuddot : float
        Pulsar rotation frequency second derivative (s^-3)

    p_orbitaldot : float
        Rate of change of orbital period

    eps1 : float
        First Laplace parameter [E x sin(lon_periest)]

    eps2 : float
        Second Laplace parameter [E x cos(lon_periest)]

    model : Binary model (BT/ELL1/DD/MSS)
    """

    def __init__(self, met0, nu0, t_orbital, p_orbital, axis_proj, ecc=0., lon_periast=0.,
                 nudot0=0., nuddot=0., p_orbitaldot=0., eps1=None, eps2=None, model=None):
        """Constructor
        """
        xEphemeris.__init__(self, met0, nu0, nudot0, nuddot)
        self.t_orbital = t_orbital
        self.p_orbital = p_orbital
        self.p_orbitaldot = p_orbitaldot
        self.ecc = ecc
        self.axis_proj = axis_proj
        self.lon_periast = lon_periast
        self.eps1 = eps1
        self.eps2 = eps2
        self.model = model

    def omega_orb(self, met):
        """Return the inverse of the orbit period at a given time.
        """
        dt = self._dt(met)
        return 1. / self.p_orbital - 0.5 * (self.p_orbitaldot / (self.p_orbital**2.)) * dt

    @classmethod
    def from_file(cls, parfile_path):
        """Read a set of attributes given a parameter file.
        """
        par_dict = read_par_file(parfile_path)
        params = {}
        params['met0'] = mjd_to_met(par_dict['PEPOCH'])
        params['nu0'] = par_dict['F0']
        params['nudot0'] = par_dict.get('F1', 0.)
        params['nuddot'] = par_dict.get('F2', 0.)
        t0 = par_dict.get('T0')
        params['t_orbital'] = mjd_to_met(par_dict.get('TASC', t0))
        params['p_orbital'] = days_to_seconds(par_dict['PB'])
        params['axis_proj'] = par_dict['A1'] # in lt-s
        params['p_orbitaldot'] = par_dict.get('PBDOT', 0.)
        params['ecc'] = par_dict.get('ECC', 0.)
        params['lon_periast'] = par_dict.get('OM', 0.)
        params['eps1'] = par_dict.get('EPS1', None)
        params['eps2'] = par_dict.get('EPS2', None)
        params['model'] = par_dict['BINARY']
        return cls(**params)

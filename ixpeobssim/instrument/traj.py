#!/urs/bin/env python
#
# Copyright (C) 2019--2020, the ixpeobssim team.
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

"""Spacecraft trajecory related functions
"""

from __future__ import print_function, division

from enum import IntEnum
import os
import numbers
import datetime

try:
    from io import BytesIO
except:
    from StringIO import StringIO as BytesIO

import numpy
import matplotlib.path
import matplotlib.patches
import skyfield
import skyfield.api
import skyfield.iokit
from skyfield.timelib import Timescale
from skyfield.iokit import parse_deltat_data
from skyfield.iokit import parse_deltat_preds, parse_leap_seconds
from skyfield.nutationlib import iau2000b
from skyfield.functions import angle_between
from skyfield.units import Angle

from ixpeobssim import IXPEOBSSIM_INSTRUMENT_DATA, IXPEOBSSIM_DATA
from ixpeobssim.utils.environment import SKYFIELD_VERSION
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.profile import timing, xChrono, xMemoryProfiler
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.utils.time_ import met_to_jd, xTimeInterval
from ixpeobssim.evt.gti import xGTIList
from ixpeobssim.instrument.mma import apply_dithering


JULIAN_DATE_UNMOD = 2400000.5
RADIUS_EARTH_AU = 4.26352e-5
AU_TO_KM = 1.495978707e8 #1 AU in km

EARTH_RADIUS = 6371.   # km
EARTH_MASS = 5.972e24  # kg
G = 6.67408e-11

_TIMESCALE_IO_ADVICE = "Try opening the same URL in your browser to learn more about the problem."


class xSkyfieldLoader(skyfield.iokit.Loader):

    def __init__(self, verbose=True, expire=True):
        """Overloaded timescale loader, pointing to IXPEOBSSIM_DATA for all the
        data files to be dowloaded locally.
        """
        skyfield.iokit.Loader.__init__(self, IXPEOBSSIM_DATA, verbose, expire)

    @staticmethod
    def _get_data(file_name):
        """Drop-in replacement function for the pkgutils.get_data() method
        used in the original timescale implementation.

        This is just picking up the file from the proper folder in the ixpeobssim
        installation, rather than from the skyfield installation.
        See https://docs.python.org/3/library/pkgutil.html for more details.
        """
        file_path = os.path.join(IXPEOBSSIM_INSTRUMENT_DATA, file_name)
        with open(file_path, 'rb') as input_file:
            data = input_file.read()
        return data

    def timescale(self, delta_t=None, builtin=False):
        """Open or download the three timescale files, returning a `Timescale`.

        This is an almost verbatim re-implementation of the timescale() method in
        skyfield/iokit.py, the only difference being the locations from where we
        retrieve the data.
        """
        logger.info('Running custom skyfield loader...')
        if builtin:
            logger.info('Using builtin timescale ancillary data (from %s)...',\
                IXPEOBSSIM_INSTRUMENT_DATA)
            b = self._get_data('deltat.data')
            # Hack to accommodate a change in skyfield 1.26.
            # See issue #309.
            if SKYFIELD_VERSION < '1.26':
                expiration_date, data = parse_deltat_data(BytesIO(b))
            else:
                data = parse_deltat_data(BytesIO(b))
            # End of hack.
            b = self._get_data('deltat.preds')
            # Hack to accommodate a change in skyfield 1.26.
            # See issue #309.
            if SKYFIELD_VERSION < '1.26':
                expiration_date, preds = parse_deltat_preds(BytesIO(b))
            else:
                preds = parse_deltat_preds(BytesIO(b))
            # End of hack.
            data_end_time = data[0, -1]
            i = numpy.searchsorted(preds[0], data_end_time, side='right')
            delta_t_recent = numpy.concatenate([data, preds[:,i:]], axis=1)
            b = self._get_data('Leap_Second.dat')
            expiration_date, arrays = parse_leap_seconds(BytesIO(b))
            leap_dates, leap_offsets = arrays

            return Timescale(delta_t_recent, leap_dates, leap_offsets)

        logger.info('Using downloaded timescale ancillary data (from %s)...',\
            self.directory)
        if delta_t is not None:
            delta_t_recent = numpy.array(((-1e99, 1e99), (delta_t, delta_t)))
        else:
            try:
                data = self('deltat.data')
                preds = self('deltat.preds')
            except IOError as e:
                e.args = (e.args[0] + _TIMESCALE_IO_ADVICE,) + e.args[1:]
                raise
            data_end_time = data[0, -1]
            i = numpy.searchsorted(preds[0], data_end_time, side='right')
            delta_t_recent = numpy.concatenate([data, preds[:,i:]], axis=1)
        try:
            leap_dates, leap_offsets = self('Leap_Second.dat')
        except IOError as e:
            e.args = (e.args[0] + _TIMESCALE_IO_ADVICE,) + e.args[1:]
            raise
        return Timescale(delta_t_recent, leap_dates, leap_offsets)



class xSAABoundary(matplotlib.path.Path):

    """Small container class representing the SAA boundary.

    The main purpose of this class is to encapsulate the algorithm to figure
    out whether any given point on the Earth surface is inside the SAA polygon
    or not. You can glance through the related issue at
    https://bitbucket.org/ixpesw/ixpeobssim/issues/232
    to get an account of the underlying discussion, but essentially we opted
    for the SAA class to be a subclass of matplotlib.path.Path since the
    latter provides all the functionality that we need (includeing plotting),
    implemented efficiently and without bringing in yet another dependence.
    """

    def __init__(self, file_name='saa_polygon.txt'):
        """Constructor.
        """
        file_path = os.path.join(IXPEOBSSIM_INSTRUMENT_DATA, file_name)
        logger.info('Loading SAA boundaries from %s...', file_path)
        vertices = numpy.loadtxt(file_path)
        vertices = numpy.vstack((vertices, vertices[0]))
        matplotlib.path.Path.__init__(self, vertices, closed=True)

    def plot(self, plot_vertices=True, **kwargs):
        """Plot the polygon.
        """
        kwargs.setdefault('facecolor', 'orange')
        kwargs.setdefault('edgecolor', 'black')
        kwargs.setdefault('alpha', 0.5)
        patch = matplotlib.patches.PathPatch(self, **kwargs)
        plt.gca().add_patch(patch)
        if plot_vertices:
            plt.plot(*numpy.hsplit(self.vertices, 2), 'o', color=kwargs.get('edgecolor'))

    def contains(self, lon, lat):
        """Return a boolean mask signaling the pairs of coordinates inside the
        polyogon.
        """
        # Support the function call with scalar arguments.
        if isinstance(lon, numbers.Number) and isinstance(lat, numbers.Number):
            return self.contains_point((lon, lat))
        # Handle array arguments---note the additional "s" in the underlying call.
        return self.contains_points(list(zip(lon, lat)))

    def __str__(self):
        """String formatting.
        """
        text = 'SAA boundary: '
        for lon, lat in self.vertices:
            text += '(%.2f, %.2f)--' % (lon, lat)
        return text.strip('-')



class xTLE:

    """Small utility function to create TLE strings for the purpose of
    our simulation.

    See https://en.wikipedia.org/wiki/Two-line_element_set for some
    background information.

    Mind this is not meant to be a complete interface, which would be silly---what
    we are after is really the ability to generate a fictional TLE for a simple
    circular orbit with arbitrary altitude and inclination. We purposedly do not
    provide the ability to write ouput file in order not to pollute the
    environment---this is something that should remain internal to the
    simulation.

    And, for completeness, here is the IXPE TLE pulled out from
    https://www.n2yo.com/satellite/?s=49954
    shortly after launch.

    1 49954U 21121A   21351.00640149  .00001120  00000-0  35770-4 0  9994
    2 49954   0.2300 281.7657 0011347 134.4260 303.9164 14.90740926  1166
    """

    @staticmethod
    def _internationl_designator(year=2021, launch_no=0, launch_piece='A'):
        """
        4 	10–11 	International Designator (last two digits of launch year)
        5 	12–14 	International Designator (launch number of the year)
        6 	15–17 	International Designator (piece of the launch)
        """
        return '%s%03d%s' % (str(year)[-2:], launch_no, launch_piece.ljust(3))


    @staticmethod
    def _epoch(year=2021, day=0.):
        """
        7 	19–20 	Epoch Year (last two digits of year)
        8 	21–32 	Epoch (day of the year and fractional portion of the day)

        For reasons that now escape me this used to be set on the basis of the
        start time of the simulation (i.e., datetime.now()), and it used to
        be different for each run, which prevented the user from being able to
        have the same GTIs twice, even if the start_date and duration of the
        observation did not change.
        https://bitbucket.org/ixpesw/ixpeobssim/issues/394

        Since the precise epoch of the TLE isn't really important---but getting
        reproducible results is---we set the epoch to a fictional point in time,
        which is the same for every run.
        """
        return '%s%012.8f' % (str(year)[-2:], day)

    @staticmethod
    def _format_inclination(inclination):
        """
        """
        if inclination <= 10.:
            return '%.4f' % inclination
        return '%.3f' % inclination

    @staticmethod
    def _mean_motion(altitude, correction_factor=0.9996733350090976):
        """Calculate the mean motion given the altitude over the earth surface.

        This is just using the elementary equation of motion for circular orbits,
        and the small correction factor was determined after the fact to fix
        a ~10 km residual error on the actual altitude.
        """
        r = 1000. * (EARTH_RADIUS + altitude)
        T = numpy.sqrt(4. * numpy.pi**2. * r**3. / G / EARTH_MASS)
        return correction_factor * 86400. / T

    @staticmethod
    def lines(inclination=0.23, altitude=601.1, catalog_no=49954, launch_year=2021,
        launch_no=121, launch_piece='A', launch_day=351.00640149, ballistic_coefficient='.00001120',
        drag_term='35770-4', element_number=9994, ascending_node=281.7657,
        eccentricity='0011347', perigee=134.4260, mean_anomaly=303.9164):
        """
        """
        intl_designator = xTLE._internationl_designator(launch_year, launch_no, launch_piece)
        epoch = xTLE._epoch(launch_year, launch_day)
        line1 = '1 %05dU %s %s  %s  00000-0  %s 0  %4d' %\
            (catalog_no, intl_designator, epoch, ballistic_coefficient, drag_term, element_number)
        inclination = xTLE._format_inclination(inclination)
        mean_motion = xTLE._mean_motion(altitude)
        line2 = '2 %05d   %s %3.4f %s %3.4f %3.4f %.8f  1166' %\
            (catalog_no, inclination, ascending_node, eccentricity, perigee, mean_anomaly, mean_motion)
        return line1, line2



class xIXPETrajectory:

    """A class to characterize the spacecraft's orbit and trajectory for
    ixpeobssim.

    This is supposed to keep track of the necessary elements to keep track
    of the state of the observatory and generate the GTI.

    The TLE data are read from the proper file in
    $IXPEOBSSIM_INSTRUMENT_DATA (which satellite precisely is controlled by
    the satellite_name argument).

    The planetary and lunar ephemeris are taken from
    `JPL <https://ssd.jpl.nasa.gov/?planet_eph_export>`_.
    While JPL distributes and maintain an extensive set of files, in order to
    avoid large downloads at runtime, we have packaged or own trimmed
    versions of DE405 and DE430t within the 2000--2040 time span into the
    ixpeobssim/instrument/data folder. The trimming has been performed using
    the jplephem utility described at
    https://rhodesmill.org/skyfield/planets.html

    The generation of the timescale information takes advantage of all the
    skyfield machinery. Note that, since the default skyfield method do not
    reasily offer any flexibility in terms of where the ancillary files are
    located, we have implemented our own Loader class that behaves as the
    the original with two notable exceptions:

    * the builtin timescale ancillary files (when necessary) are loaded from
      IXPEOBSSIM_INSTRUMENT_DATA, rather than from the sjyfield installation;
    * the updated timescale ancillary files (when necessary) are downloaded
      (and read from) IXPEOBSSIM_DATA.

    Args
    ----
    ephemeris_file_name : str
        The ephemeris file name (in the IXPEOBSSIM_INSTRUMENT_DATA folder).

    builtin_timescale: bool
        Flag passed to the skyfield timescale generation routine (if True the
        static files shipped with ixpeobssim are used, otehrwise they are
        downloaded).
    """

    def __init__(self, ephemeris_file_name='de430t_2000_2040.bsp',
                 min_sun_angle=65., max_sun_angle=115., builtin_timescale=True):
        """Constructor.
        """
        self.min_sun_angle = min_sun_angle
        self.max_sun_angle = max_sun_angle
        # Load the timescale.
        self.timescale = xSkyfieldLoader().timescale(builtin=builtin_timescale)
        # Load the SAA boundary.
        self.saa_boundary = xSAABoundary()
        # Load the satellite TLE.
        self.satellite = skyfield.api.EarthSatellite(*xTLE.lines(), 'IXPE', self.timescale)
        logger.debug(self.satellite)
        # Load the planetary ephemeris.
        ephem_file_path = os.path.join(IXPEOBSSIM_INSTRUMENT_DATA, ephemeris_file_name)
        logger.info('Loading planet ephemeris from %s...', ephem_file_path)
        self.planets = skyfield.api.load_file(ephem_file_path)

    @staticmethod
    def load_tle(tle_file_name='tle_data_science.txt', satellite_name='AGILE'):
        """
        Load TLE data for a given satellite.
        Note the tle() method is deprecated with an entertaining comment
        DEPRECATED: in a misguided attempt to be overly convenient, this routine
        builds an unweildy dictionary of satellites with keys of two different
        Python types: integer keys for satellite numbers, and string keys for
        satellite names. It even lists satellites like ``ISS (ZARYA)`` twice, in
        case the user wants to look them up by a single name like ``ZARYA``.
        What a mess. Users should instead call the simple ``tle_file()`` method,
        and themselves build any dictionaries they need.

        Args
        ----
        tle_file_name : str
            The name of the file (in the IXPEOBSSIM_INSTRUMENT_DATA folder)
            containing the TLE data.

        satellite_name : str
            The name of the satellite to be looked for in the TLE data ancillary
            files. In absence of sensible TLE data for IXPE, NUSTAR and AGILE
            are two scientific satellites in near-equatorial orbit that can be
            used for our purpose.
        """
        tle_file_path = os.path.join(IXPEOBSSIM_INSTRUMENT_DATA, tle_file_name)
        logger.info('Loading IXPE TLE data from %s...', tle_file_path)
        satellites = {s.name: s for s in skyfield.api.load.tle_file(tle_file_path)}
        try:
            return satellites[satellite_name]
        except KeyError:
            abort('Cannot retrieve TLE data for satellite %s' % satellite_name)

    def __del__(self):
        """This is to close the underlying ephemeris file---I am even unsure
        this is a good idea.

        See https://github.com/skyfielders/python-skyfield/issues/374
        for more information.
        """
        self.planets.spk.close()

    def met_to_ut(self, met, use_iau2000b=True):
        """Convert MET(s) to skyfield.Time object(s), using the proper
        underlying timescale.

        See https://github.com/skyfielders/python-skyfield/issues/373
        for all the gory details of the iau2000b magic! That really saved the day.
        """
        ut = self.timescale.ut1_jd(met_to_jd(met))
        if use_iau2000b:
            ut._nutation_angles = iau2000b(ut.tt)
        return ut

    def geocentric(self, met):
        """Return the geocentric position at a given MET.
        """
        ut = self.met_to_ut(met)
        return self.satellite.at(ut)

    def position(self, met):
        """Return the position (longitude and latitude in decimal degrees) of
        the satellite at a given MET.
        """
        pos = self.geocentric(met).subpoint()
        return pos.longitude.degrees, pos.latitude.degrees

    def elevation(self, met):
        """Return the altitude of the satellite (in km) at a given MET.
        """
        return self.geocentric(met).subpoint().elevation.km

    def in_saa(self, met):
        """Return whether the satellite is in the SAA at a given time.
        """
        lon, lat = self.position(met)
        return self.saa_boundary.contains(lon, lat)

    def target_planet_angle(self, met, target_ra, target_dec, planet_name):
        """Return the angular separation between a target RA and DEC and a
        planet at a given MET time.

        Warning
        -------
        For now the angular separation is calculated from a geocentric (i.e.,
        not spacecraft) position, which should be sufficient for our needs.
        """
        # Support the call for scalar met values.
        if isinstance(met, numbers.Number):
            met = numpy.array([met])
        ut = self.met_to_ut(met)
        # RA for target_pos needs to be given in RA hours, not degrees
        # (15 degrees per hour of RA)
        # target_pos = skyfield.api.position_of_radec(target_ra / 15., target_dec)
        # The above version is for a future version of skyfield with better
        # documentation and a more sensible default distance.
        # It is implemented identically using the current (as of 1.20) version below.
        # The distance argument corresponds to 1 Gpc in AU.
        target_pos = skyfield.api.position_from_radec(target_ra / 15., target_dec,
                                                      distance=206264806247096.38)
        planet_pos = self.planets['earth'].at(ut).observe(self.planets[planet_name])

        # Hack to maintain backward compatibility with older versions of skyfield.
        if SKYFIELD_VERSION < '1.23':
            return planet_pos.separation_from(target_pos).degrees
        # In newer version of skyfield we seem to have broadcasting issues, see
        # https://bitbucket.org/ixpesw/ixpeobssim/issues
        target_pos = target_pos.position.au
        planet_pos = planet_pos.position.au
        if target_pos.ndim == 1 and planet_pos.ndim == 2:
            target_pos = target_pos.reshape((3, 1))
        return Angle(radians=angle_between(target_pos, planet_pos)).degrees

    def target_sun_angle(self, met, target_ra, target_dec):
        """Return the angular separation between a target RA and DEC and the Sun
        at a given MET.
        """
        return self.target_planet_angle(met, target_ra, target_dec, 'sun')

    def target_moon_angle(self, met, target_ra, target_dec):
        """Return the angular separation between a target RA and DEC and the
        Moon at a given MET.
        """
        return self.target_planet_angle(met, target_ra, target_dec, 'moon')

    def sun_constrained(self, met, target_ra, target_dec):
        """Return whether a given target at a given time cannot be observed
        due to the Sun constraints.
        """
        sun_angle = self.target_sun_angle(met, target_ra, target_dec)
        return numpy.logical_or(sun_angle < self.min_sun_angle, sun_angle > self.max_sun_angle)

    def target_occulted(self, met, target_ra, target_dec, limiting_altitude = 200.):
        """Return whether a given target is occulted at a given time.
        """
        # Support the call for scalar met values.
        if isinstance(met, numbers.Number):
            met = numpy.array([met])

        ut = self.met_to_ut(met)
        #This target position is by default at a distance of 1 gigaparsec
        #Much too far to worry about parallax or other effects but also far enough
        #to cause floating point precision errors if the vector algebra isn't done carefully
        #target_pos = skyfield.positionlib.position_of_radec(target_ra / 15., target_dec,  t=ut,epoch=self.timescale.J2000)
        #The distance corresponds to 1 Gpc in AU
        target_pos = skyfield.api.position_from_radec(target_ra / 15., target_dec, distance = 206264806247096.38  )

        target_ecliptic = target_pos.ecliptic_position().au


        #The target_ecliptic position has a total length/distance of 1e6 AU.
        #which corresponds to ~ 4.8 kpc
        #The actual distance to the source past 4.8 pc is irrelevant to
        #the occultation calculation

        earth_ecliptic = self.planets['earth'].at(ut).ecliptic_position().au
        sat_ecliptic = self.geocentric(met).ecliptic_position().au


        limit_dist = (EARTH_RADIUS + limiting_altitude) / AU_TO_KM
        #print(limit_dist)
        dist_perp_list = []
        dist_sattarg_list = []
        dist_targearth_list = []
        ttt_list = []
        n_times = numpy.shape(ut)[0]
        #At each time we calculate the minimum perpendicular distance
        #from the Earth position to the vector connecting
        #target ra/dec and satellite
        #Parameterizing the vector
        #v = sat + elciptic * t
        #where t is an arbitary parameter
        #we can calculate the minimum perpendicular distance to the Earth position
        #and value of t when it reaches that position (denoted at ttt in the code).
        #If the minimum perpendicular distance is within R_Earth + limit_altitude
        #and ttt > 0 (meaning the minimum distance is reached as a vector goes from
        #satellite to target then the Earth is occulting the observation.
        #All math is given at this website
        #https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        for j in range(n_times):

            this_earth = earth_ecliptic[:,j]
            this_sat = sat_ecliptic[:,j] + this_earth
            # sat_ecliptic is natively a geocentric coordinate, while the others
            #are are Solar Sys baycentric


            diff_earth_sat = this_earth - this_sat
            l_earth_sat = numpy.sqrt(numpy.dot(diff_earth_sat,diff_earth_sat))

            diff_earth_targ = this_earth - target_ecliptic
            l_earth_targ = numpy.sqrt(numpy.dot(diff_earth_targ,diff_earth_targ))

            diff_targ_sat = target_ecliptic - this_sat
            l_targ_sat = numpy.sqrt(numpy.dot(diff_targ_sat, diff_targ_sat))


            #This is the (x0-x1) X ( x0-x2) version of the min dist formula
            cross = numpy.cross(diff_earth_sat,diff_earth_targ)

            lcross = numpy.sqrt(numpy.dot(cross,cross))

            dist_perp = lcross / l_targ_sat

            #Floating point precision is an issue for this calculation since the
            # diff_eath_sat << diff_earth_targ. So we actually calculate the value
            # of the parameter t where it reaches minimum distance.

            ttt = -1.* numpy.dot( -1. * diff_earth_sat, diff_targ_sat) / l_targ_sat**2


            dist_perp_list.append(dist_perp)

            ttt_list.append(ttt)


        dist_perp_arr = numpy.array(dist_perp_list)
        ttt_arr = numpy.array(ttt_list)

        #Minimum perpendicular distance between Earth position and target RA/DEC
        #crosses the Earth
        perp_mask = dist_perp_arr < limit_dist

        #The Earth is closer to the target RA/DEC than the satellite
        front_mask = ttt_arr > 0.

        #True means is occulted, and False means the RA/DEC is not occulted.
        return numpy.logical_and(perp_mask,front_mask)

    @staticmethod
    def equispaced_met_grid(start_met, stop_met, approximate_step):
        """Return an equally-spaced grid of MET values according to the input
        argument.

        Note that the array returned ranges exactly from start_met to stop_met,
        and the step size is rounded to guarantee that.
        """
        num_steps = int((stop_met - start_met) / approximate_step) + 1
        return numpy.linspace(start_met, stop_met, num_steps)

    @staticmethod
    def _interleave_arrays(a1, a2):
        """Interleave two numpy arrays, inserting the midpoint for each pair.

        This is used for our poor-man binary search of the SAA transit times and
        the appearance/disappearance of the target source. We typically pass two
        series of MET values, and return the interleaved array, as well as the
        maximum "error" on the midpoint, i.e., the absolute difference between
        the highest and lowest value in each input pair.

        Shamelessly inspired from https://stackoverflow.com/questions/5347065
        """
        assert a1.size == a2.size
        assert (a2 > a1).all()
        out = numpy.empty((3 * a1.size,), dtype=a1.dtype)
        out[0::3] = a1
        out[1::3] = 0.5 * (a1 + a2)
        out[2::3] = a2
        out = numpy.unique(out)
        max_error = (0.5 * (a2 - a1)).max()
        return out, max_error

    def _generic_binary_search(self, start_met, stop_met, function, args=(),
                               initial_step=100., precision=0.001):
        """Generic binary search to find transitions in time into a function of
        the satellite coordinates.

        Args
        ----
        met : array_like
            The initial array of MET values.

        pos : 2-element tuple (lon, lat) of array_like
            The initial tuple of position arrays.

        function : callable
            The function(lon, lat) we want to apply the binary search to.

        args : tuple
            The optional additional arguments to be passed to the function.

        precision : float
            The minimum precision (in s) to which we want to calculate the
            transition times.
        """
        met = self.equispaced_met_grid(start_met, stop_met, initial_step)
        err = 2 * precision
        while err >= precision:
            # Apply the releant function for the binary search.
            mask = function(met, *args)
            # Calculate the indices where the value of the function flips from
            # True to False or vice-versa.
            idx = numpy.argwhere(numpy.diff(mask)).squeeze()
            # Handle the case where the index mask is empty---i.e. the are no
            # flips in the objective function over the met array. This means that
            # we either return the entire MET interval or nothing.
            if idx.size == 0:
                # The function is true at the beginning and true for the entire
                # [start_met, stop_met] interval.
                if mask[0] == True:
                    return [(start_met, stop_met)]
                # The function is never True in the parent interval.
                return []
            # Create a new array with only the MET values around the flips, adding
            # the midpoint for the next step of the binary search.
            # At this point met is a one-dimensional array containing triplets
            # (low, mid, high) of values, where mid is our best estimate of the
            # transition and err is the maximum of (high - low), representing
            # the worst case maximum error.
            met, err = self._interleave_arrays(met[idx], met[idx + 1])
        # Take the midpoints of the final triplets--that's all we care about.
        met = met[1::3]
        # We now evaluate the position and the mask *right after* the transition
        # to figure out whether we are entering or exiting.
        entrance = function(met + precision, *args)
        # Take care of the corner cases where the first transition is an exit or
        # the last transition is an entrance. Note that we cannot test the
        # equality of the elements of a numpy boolean array with is, need to use ==.
        if entrance[0] == False:
            met = numpy.insert(met, 0, start_met)
        if entrance[-1] == True:
            met = numpy.append(met, stop_met)
        return met

    @staticmethod
    def _array_to_epochs(met):
        """Small convenince function to transform an array of MET values into
        a list of epoch, i.e., a list of (start, stop) two-elements tuples.

        This is turning, e.g., [1., 2., 3.] into [(1., 2.), (2., 3.), (3., 4.)]
        """
        return [(met[i], met[i + 1]) for i in range(0, len(met), 2)]

    def saa_mets(self, start_met, stop_met, initial_step=100., precision=0.001):
        """Return the list of METs for which we cross the SAA boundaries within
        a given time interval.
        """
        args = start_met, stop_met, self.in_saa, (), initial_step, precision
        return self._generic_binary_search(*args)

    def saa_epochs(self, start_met, stop_met, initial_step=100., precision=0.001):
        """Calculate the SAA epochs.
        """
        args = start_met, stop_met
        return self._array_to_epochs(self.saa_mets(*args, initial_step, precision))

    def earth_occultation_mets(self, start_met, stop_met, target_ra, target_dec,
                               initial_step=100., precision=0.001):
        """Return the list of METS for which the Earth occultation status is
        changing.
        """
        args = start_met, stop_met, self.target_occulted, (target_ra, target_dec)
        return self._generic_binary_search(*args, initial_step, precision)

    def earth_occultation_epochs(self, start_met, stop_met, target_ra, target_dec,
                                 initial_step=100., precision=0.001):
        """Calculate the epochs when a given target is occulted by the Earth.
        """
        args = start_met, stop_met, target_ra, target_dec, initial_step, precision
        return self._array_to_epochs(self.earth_occultation_mets(*args))

    def sun_constraint_mets(self, start_met, stop_met, target_ra, target_dec,
                            initial_step=1000., precision=0.001):
        """Calculate the list of METs for which the Sun observability changes.
        """
        args = start_met, stop_met, self.sun_constrained, (target_ra, target_dec)
        return self._generic_binary_search(*args, initial_step, precision)

    def sun_constraint_epochs(self, start_met, stop_met, target_ra, target_dec,
                              initial_step=1000., precision=0.001):
        """Calculate the epochs when observations are inhibited by the Sun
        constraints.
        """
        args = start_met, stop_met, target_ra, target_dec, initial_step, precision
        return self._array_to_epochs(self.sun_constraint_mets(*args))

    def timeline_mets(self, start_met, stop_met, target_ra, target_dec, saa, occult,
                      initial_step=100., precision=0.001):
        """Return a list of all the relevant MET values for the calculation of
        the observation timeline, i.e., those that pertain to the SAA and to
        the Earth occultation.

        The calculation is started with the observaton start and stop met, and
        the SAA and Earth occultation MET values are addedd if necessary.

        The function returns separate lists for the SAA and Earth occultation
        MET values, as well as a sorted logical OR of the two with no, diplicates,
        since all is needed downstream to calculate the timeline.
        """
        mets = numpy.array([start_met, stop_met])
        saa_mets = numpy.array([])
        occult_mets = numpy.array([])
        if saa:
            args = start_met, stop_met, initial_step, precision
            saa_mets = self.saa_mets(*args)
            mets = numpy.append(mets, saa_mets)
        if occult:
            args = start_met, stop_met, target_ra, target_dec, initial_step, precision
            occult_mets = self.earth_occultation_mets(*args)
            mets = numpy.append(mets, occult_mets)
        return numpy.unique(mets), saa_mets, occult_mets

    @staticmethod
    def complement_epochs(epochs, start_met, stop_met):
        """Return the complement to a list of epochs.
        """
        # If the list of epochs is empty, then just return a list with a single
        # epoch covering the entire observation.
        if len(epochs) == 0:
            return [(start_met, stop_met)]
        met = numpy.array(tuple(_met for epoch in epochs for _met in epoch))
        if met[0] > start_met:
            met = numpy.insert(met, 0, start_met)
        else:
            met = met[1:]
        if met[-1] < stop_met:
            met = numpy.append(met, stop_met)
        else:
            met = met[:-1]
        return [(met[i], met[i + 1]) for i in range(0, len(met), 2)]

    @staticmethod
    def _epochs_to_arrays(epochs):
        """Tranform a list of epochs, i.e., a list of length-2 tuples of the
        form (t1, t2) into a pair of arrays suitable for the determination of
        the GTIs.

        .. warning:
            This is only used in the gti_list() class method, which is now
            deprecated, and can be removed when, and if, the latter is removed.
        """
        # The first array contains all the epoch bounds, unpacked from the
        # original tuples.
        met = numpy.array(tuple(_met for epoch in epochs for _met in epoch))
        # The second array signal whether each of the MET value signal an
        # enter into the epoch (+1) or an exit from the epoch (-1).
        sign = numpy.array([1, -1] * len(epochs))
        return met, sign

    def gti_list(self, start_met, stop_met, target_ra, target_dec,
                 saa=True, occult=True, initial_step=100., precision=0.001):
        """Calculate the GTIs, i.e., the ephocs where we are not in the SAA
        and the target is not occulted by the Earth.

        While this method seems (and is, in fact) quite convoluted, we decided
        to calculate separately the SAA and Earth occultation epochs, when the
        instrument is not taking data, and merge them at the end to extract the
        good time intervals. The reasoning is that, while the former epochs are
        quite regular over the timescale of an observation---which allows to
        perform a sensible binary search starting from a reasonably corse grid
        (say with a 100 s step)---when we put in AND the two we get all kind of
        resonances, to the point that the GTIs can be arbitrarily small, and
        we would need a much finer (and computationally intensive) initial grid
        in order not to introduce measurable errors.

        .. warning:
            The xObservationTimeline class provide a fresh rewite of the GTI
            logics, in a simpler way that allows to keep track of the precise
            status of the observatory at each point in time. This function
            is temporarely kept for debugging reason, but will be effectively
            removed.

        Args
        ----
        start_met : float
            The initial MET for the observation.

        stop_met : float
            The final MET for the observation.

        target_ra : float
            The right ascension of the target source in decimal degrees.

        target_dec : float
            The declination of the target source in decimal degrees.

        saa : bool
            Consider the SAA passages in the GTI calculation.

        occult : bool
            Consider the Earth occultation in the GTI calculation.

        initial_step : float
            The step (in s) for the initial equispacedtime time grid used to
            seed the binary search. (This can slow things down when too small.)

        precision : float
            The worst-case precision of the GTI calculation (i.e., the stop
            condition for the binary search)
        """
        # Trivial case: we return a list containing a single GTI covering the
        # entire observation :-)
        saa_epochs = []
        occult_epochs = []
        if saa is False and occult is False:
            return xGTIList(start_met, stop_met, (start_met, stop_met)), saa_epochs, occult_epochs
        # If we're here we got something to do. Prepare the data structures.
        met = numpy.array([], dtype=float)
        sign = numpy.array([], dtype=int)
        epochs_list = []
        # Collect the information we need about the parts of the observation
        # to be trimmed out.
        if saa is True:
            args = start_met, stop_met, initial_step, precision
            saa_epochs = self.saa_epochs(*args)
            epochs_list.append(saa_epochs)
        if occult is True:
            args = start_met, stop_met, target_ra, target_dec, initial_step, precision
            occult_epochs = self.earth_occultation_epochs(*args)
            epochs_list.append(occult_epochs)
        # Loop over the list of epochs where we are not taking data and
        # merge them into a unique array of MET values, keeping track of the
        # enter/exit flags.
        for epochs in epochs_list:
            _met, _sign = self._epochs_to_arrays(epochs)
            met = numpy.append(met, _met)
            sign = numpy.append(sign, _sign)
        # Handle the case where the MET list is empty---this means that we
        # have a single time interval covering the entire observation.
        if met.size == 0:
            return xGTIList(start_met, stop_met, (start_met, stop_met)), saa_epochs, occult_epochs
        # Sort the MET values and the signs
        idx = numpy.argsort(met)
        met = met[idx]
        sign = sign[idx]
        # We still have the potential issue of overlapping MET values, e.g.,
        # when at the end of the observation we are occulted *and* in the SAA,
        # so that stop_met occurs twice at the end of the array.
        met, idx = numpy.unique(met, return_index=True)
        sign = sign[idx]
        # Final loop over the transitions. We keep track of the sum of the sign
        # values, and every time we hit zero, we're entering a GTI.
        sign_sum = 0
        gti_list = xGTIList(start_met, stop_met)
        # If the first epoch entrance is after the start MET, we have a GTI
        # right away.
        if met[0] > start_met:
            gti_list.append_gti(start_met, met[0])
        # If the last epoch exit is before the stop MET, then we need to append
        # the stop MET itself to the MET list in order to avoid an index error
        # at the last step of the loop.
        if met[-1] < stop_met:
            met = numpy.append(met, stop_met)
        for i, _sign in enumerate(sign):
            sign_sum += _sign
            if sign_sum == 0 and i < len(met) - 1:
                gti_list.append_gti(met[i], met[i + 1])
        return gti_list, saa_epochs, occult_epochs



class xTimelineEpoch(xTimeInterval):

    """Small convenience class to encapsulate an "epoch" within a given
    observation.

    In this context an epoch is the datum of a start_met, a stop_met, and a
    series of bit flags indicating whether the observatory is, e.g., in the SAA
    or occulted by the Earth.

    See https://bitbucket.org/ixpesw/ixpeobssim/issues/417
    """

    def __init__(self, start_met, stop_met, in_saa, occulted):
        """Constructor.
        """
        xTimeInterval.__init__(self, start_met, stop_met)
        self.in_saa = in_saa
        self.occulted = occulted

    def bitmask(self):
        """Return a bit mask corresponding to the observatory status during the
        epoch.
        """
        return self.in_saa << 1 | self.occulted

    def shrink(self, start_padding, stop_padding):
        """Create a new epoch where the bounds are shrinked by given amounts on
        either side. The "shrinked" epoch is inhering the SAA and occultation flags.

        Note that the original object is left untouched, i.e., the new epoch is
        created in place.
        """
        start_met = self.start_met + start_padding
        stop_met = self.stop_met - stop_padding
        return self.__class__(start_met, stop_met, self.in_saa, self.occulted)

    def isgti(self):
        """Return True if the epoch is a good time interval, i.e., we're not in
        the SAA and we're not occulted by the Earth.
        """
        return not self.in_saa and not self.occulted

    def isocti(self):
        """Return True of the epoch is a good on-orbit calibration time interval,
        i.e., we are occulted by the Earth but not in the SAA.
        """
        return self.occulted and not self.in_saa

    def __str__(self):
        """String formatting.
        """
        return 'Epoch %s [b{:02b}]'.format(xTimeInterval.__str__(self), self.bitmask())



class xObservationTimeline(xTimeInterval):

    """Small container class encapsulating the detailed planning, in terms of
    GTIs, calibration intervals and SAA passages, for a specific observation
    (i.e., a portion of the mission with a specific, fixed pointing).

    This is complete rewrite of the code that used to live into the gti_list()
    method of the xIXPETrajectory class, ini response to
    https://bitbucket.org/ixpesw/ixpeobssim/issues/417
    The xIXPETrajectory class is still responsible for all of the complex logics
    and calculation to detect the relevant MET values for the SAA passages and
    the Earth occultation, and these values are turned into a list of
    xTimelineEpoch objects to be consumed downstream.

    This class provides a coherent and unified interface to extract good time
    intervals and onboard calibration intervals, with optional padding on
    either end and with the possibility of enforcing a minimum duration
    (after the padding).
    """

    def __init__(self, start_met, stop_met, target_ra, target_dec, saa, occult):
        """Constructor.
        """
        xTimeInterval.__init__(self, start_met, stop_met)
        self.target_ra = target_ra
        self.target_dec = target_dec
        self.trajectory = xIXPETrajectory()
        args = start_met, stop_met, target_ra, target_dec, saa, occult
        self.epochs = self._calculate_epochs(*self.trajectory.timeline_mets(*args))

    def epoch_data(self):
        """Return the epoch data in a form that can be used to fill the TIMELINE
        extension in the output FITS files.
        """
        start = numpy.array([epoch.start_met for epoch in self.epochs])
        stop = numpy.array([epoch.stop_met for epoch in self.epochs])
        in_saa = numpy.array([epoch.in_saa for epoch in self.epochs])
        occulted = numpy.array([epoch.occulted for epoch in self.epochs])
        return (start, stop, in_saa, occulted)

    @staticmethod
    def _align_met(met, grid_step, ceil=False):
        """Align a MET value, or an array of MET value, with a pre-existing
        grid of fixed steps.

        This is used to generate the SC data on a regular time gird with a
        given step.
        """
        met = int(met / grid_step) * grid_step
        if ceil:
            met += grid_step
        return met

    def sc_data(self, time_step, dither_params=None, saa=True, occult=True):
        """Return the spacecraft data to be written into the SC_DATA optional
        extension of the output photon list.

        Spacecraft data are instantaneous data (e.g., longitude, latitude,
        elevation) evaluated on a regular grid with constant spacing. Rather than
        with the start or the stop of the observation, the grid is aligned
        with MET = 0, and padded as necessary on both ends to include the
        entire span of the observation.

        The relevant data are returned into the form of a Python dictionary of the
        form {column_name: value_array} that can be used directly to fill the
        columns of the output binary table.


        Arguments:
        -----------
        time_step: float
            Interval used to save the SC information

        dither_params : (amplitude, pa, px, py) tuple, optional
            The parameters for the dithering of the observatory. The dithering is
            not applied if this is set to None.

        saa : bool
            Flag indicating whether the time intervals in the SAA should be
            flagged---the corresponding column with be identically zero if this
            is false.

        occult : bool
            Flag indicating whether the time intervals when the target is occulted
            should be flagged---the corresponding column with be identically zero
            if this is false.
        """
        met0 = self._align_met(self.start_met, time_step)
        met1 = self._align_met(self.stop_met, time_step, ceil=True)
        met_grid = numpy.arange(met0, met1 + 0.5 * time_step, time_step)
        ra_pnt, dec_pnt = apply_dithering(met_grid, self.target_ra, self.target_dec, dither_params)
        lon, lat = self.trajectory.position(met_grid)
        elevation = self.trajectory.elevation(met_grid)
        sun_angle = self.trajectory.target_sun_angle(met_grid, ra_pnt, dec_pnt)
        if saa:
            in_saa = self.trajectory.in_saa(met_grid)
        else:
            in_saa = numpy.zeros(met_grid.shape, dtype=int)
        if occult:
            args = met_grid, self.target_ra, self.target_dec
            target_occulted = self.trajectory.target_occulted(*args)
        else:
            target_occulted = numpy.zeros(met_grid.shape, dtype=int)
        return {'MET': met_grid, 'RA_PNT': ra_pnt, 'DEC_PNT': dec_pnt,
            'LAT_GEO': lat, 'LON_GEO': lon, 'ALT_GEO': elevation, 'SUN_ANGLE': sun_angle,
            'IN_SAA': in_saa, 'TARGET_OCCULT': target_occulted}

    @staticmethod
    def _bisect_odd(array_, value):
        """Return whether a single value is bisecting a given array in a odd index.

        This is useful to figure out whether we are in the SAA and/or occulted at
        any given time, by just bisecting the corresponding array of MET values.
        """
        return numpy.searchsorted(array_, value) % 2 == 1

    @staticmethod
    def _calculate_epochs(mets, saa_mets, occult_mets):
        """Loop over the relevant MET values and populate the list of epochs
        for the observation.

        Note that, for each ecpoch, we assign the SAA and Earth occultation flags
        by bisecting the corresponding MET list with the epoch center.
        """
        epochs = []
        for i, start in enumerate(mets[:-1]):
            stop = mets[i + 1]
            met0 = 0.5 * (start + stop)
            in_saa = xObservationTimeline._bisect_odd(saa_mets, met0)
            occulted = xObservationTimeline._bisect_odd(occult_mets, met0)
            epochs.append(xTimelineEpoch(start, stop, in_saa, occulted))
        return epochs

    def filter_epochs(self, min_duration=0., start_padding=0., stop_padding=0.):
        """Filter the epochs by duration.
        """
        assert min_duration >= 0.
        min_duration += start_padding + stop_padding
        return [epoch for epoch in self.epochs if epoch.duration > min_duration]

    def gti_list(self, min_duration=0., start_padding=0., stop_padding=0.):
        """Return the list of the good time intervals (GTI).

        The GTIs are defined as those when we are not in the SAA and we are
        not occulted by the Earth.

        Note the return value of this function is not a plain old list, but a
        xGTIList object.
        """
        epochs = self.filter_epochs(min_duration, start_padding, stop_padding)
        gtis = [epoch.bounds() for epoch in epochs if epoch.isgti()]
        return xGTIList(self.start_met, self.stop_met, *gtis)

    def octi_list(self, min_duration=0., start_padding=0., stop_padding=0.):
        """Return the list of the onboard calibration time intervals (OCTI).

        The OCTIs are defined as those when we are occulted by the Earth, but
        outside the SAA.
        """
        epochs = self.filter_epochs(min_duration, start_padding, stop_padding)
        return [epoch.bounds() for epoch in epochs if epoch.isocti()]

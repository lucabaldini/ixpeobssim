#!/urs/bin/env python
#
# Copyright (C) 2018, the ixpeobssim team.
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

"""Module for time conversions and related utilities.
"""

from __future__ import print_function, division


import time
import datetime

from matplotlib.dates import date2num


# pylint: disable=invalid-name




# The mission start date and time.
MISSION_START_DATETIME = datetime.datetime(2017, 1, 1)

# The Unix time of the mission start (January 1, 2017).
# This is shamelessly taken from http://www.unixtimestamp.com
MISSION_START_UNIX_TIME = 1483228800

# Modified Julian Date at the mission start (January 1, 2017).
# Taken from http://www.csgnetwork.com/julianmodifdateconv.html
MISSION_START_MJD = 57754

# Fractional part of the reference MJD (32.184 secs + 37 leap secs).
MISSION_START_MJDREFF = 8.0074074074e-04

# Default datetime format string.
DATETIME_FMT = '%Y-%m-%dT%H:%M:%S.%f'

# Launch date.
LAUNCH_DATETIME = '2021-12-09T06:00:00.0'



class xTimeInterval:


    """Small convenience class to encapsulate the very concept of a time interval.

    This was added after https://bitbucket.org/ixpesw/ixpeobssim/issues/417
    in an attempt to avoid code duplications wherever we have objects (e.g., GTIs,
    observation epochs or calibration runs) that have a start and a stop time.

    Note that all the times are assumed to be in MET.
    """

    def __init__(self, start_met, stop_met):
        """Constructor.
        """
        if stop_met <= start_met:
            raise ValueError
        self.start_met = start_met
        self.stop_met = stop_met

    def bounds(self):
        """Return the bounds of the epoch in the form of a two-element tuple (start, stop).
        """
        return (self.start_met, self.stop_met)

    @property
    def duration(self):
        """Return the total duration of the time interval in seconds.
        """
        return self.stop_met - self.start_met

    def __str__(self):
        """String formatting.
        """
        return '%.3f--%.3f (%.3f s)' % (self.start_met, self.stop_met, self.duration)



class UTCTimezone(datetime.tzinfo):

    """Derived tzinfo object to support the UTC timezone with Python 2.

    See https://docs.python.org/2/library/datetime.html#tzinfo-objects
    for more details.
    """

    __ZERO = datetime.timedelta(0)

    def utcoffset(self, dt):
        return self.__ZERO

    def dst(self, dt):
        return self.__ZERO

    def tzname(self, dt):
        return 'UTC'


"""UTC tzinfo object.
"""
try:
    # In Python 3 we have native support for UTC tzinfo.
    UTC_TZINFO = datetime.timezone.utc
except AttributeError:
    # Hack for Python 2: use the custom tzinfo derived class defined above.
    UTC_TZINFO = UTCTimezone()


def unix_to_met(ut):
    """Convert a Unix time to a MET.

    Args
    ----
    ut : float
        The input Unix time.

    Returns
    -------
    float
        The mission elapsed time corresponding to the input Unix time.
    """
    return ut - MISSION_START_UNIX_TIME


def met_to_unix(met):
    """Convert a MET to a Unix time.

    Args
    ----
    met : float
        The input mission elapsed time.

    Returns
    -------
    float
        The Unix time corresponding to the input mission elapsed time.
    """
    return met + MISSION_START_UNIX_TIME


def _datetime(ut, tzinfo=UTC_TZINFO):
    """Convenience "private" function to convert a Unix time into a
    datetime object.

    Args
    ----
    ut : float
        The input Unix time.

    tzinfo : a datetime.timezone instance or None
        The timezone info (use None for local time).

    Returns
    -------
    datetime.datetime instance
        The date and time corresponding to the input Unix time.
    """
    return datetime.datetime.fromtimestamp(ut, tzinfo)


def unix_to_string(ut, tzinfo=UTC_TZINFO, fmt=DATETIME_FMT):
    """Convert a Unix time to a string expressing time and date.

    Args
    ----
    ut : float
        The input Unix time.

    tzinfo : a datetime.timezone instance or None
        The timezone info (use None for local time).

    fmt : string
        The format for the output string

    Returns
    -------
    string
        The string corresponding to the input time
    """
    return _datetime(ut, tzinfo).strftime(fmt)


def unix_to_string_utc(ut, fmt=DATETIME_FMT):
    """Convert a Unix time to a string representing UTC date and time.

    Args
    ----
    ut : float
        The input Unix time.

    fmt : string
        The format for the output string

    Returns
    -------
    string
        The string corresponding to the input time
    """
    return unix_to_string(ut, UTC_TZINFO, fmt)


def unix_to_string_local(ut, fmt=DATETIME_FMT):
    """Convert a Unix time to a string representing local date and time.

    Args
    ----
    ut : float
        The input Unix time.

    fmt : string
        The format for the output string

    Returns
    -------
    string
        The string corresponding to the input time
    """
    return unix_to_string(ut, None, fmt)


def met_to_string(met, tzinfo=UTC_TZINFO, fmt=DATETIME_FMT):
    """Convert a MET to a string expressing time and date.

    Args
    ----
    met : float
        The input mission elapsed time.

    tzinfo : a datetime.timezone instance or None
        The timezone info (use None for local time).

    fmt : string
        The format for the output string

    Returns
    -------
    string
        The string corresponding to the input time
    """
    return unix_to_string(met_to_unix(met), tzinfo, fmt)


def met_to_string_utc(met, fmt=DATETIME_FMT):
    """Convert a MET to a string representing UTC date and time.

    Args
    ----
    met : float
        The input mission elapsed time.

    fmt : string
        The format for the output string

    Returns
    -------
    string
        The string corresponding to the input time
    """
    return met_to_string(met, UTC_TZINFO, fmt)


def met_to_string_local(met, fmt=DATETIME_FMT):
    """Convert a MET to a string representing local date and time.

    Args
    ----
    met : float
        The input mission elapsed time.

    fmt : string
        The format for the output string

    Returns
    -------
    string
        The string corresponding to the input time
    """
    return met_to_string(met, None, fmt)


def _datetime_to_timestamp(dt, tzinfo):
    """Python 2 emulation of datetime.datetime.timestamp() method.

    This is used with Python 2x in the _timestamp() convenience function below.
    See https://stackoverflow.com/questions/19801727 and
    https://hg.python.org/cpython/file/3.3/Lib/datetime.py#l1428 for
    more information.
    """
    if tzinfo is None:
        return time.mktime(dt.timetuple()) + dt.microsecond / 1.e6
    _epoch = datetime.datetime(1970, 1, 1, tzinfo=tzinfo)
    return (dt - _epoch).total_seconds()


def _timestamp(string, tzinfo=UTC_TZINFO, fmt=DATETIME_FMT):
    """Convenience "private" function to convert a string representing a
    date and time into a Unix time.

    Args
    ----
    string : string
        The input datetime string.

    tzinfo : a datetime.timezone instance or None
        The timezone info (use None for local time).

    fmt : string
        An optional format specifier for non standard input string

    Returns
    -------
    float
        The Unix time corresponding to the input string.
    """
    # Create a naive datetime object corresponding to the input string.
    dt = datetime.datetime.strptime(string, fmt)
    # If necessary set the timezone information (note that datetime objects
    # are immutable, so we need to create a new copy, here.)
    if dt is not None:
        dt = dt.replace(tzinfo=tzinfo)
    # And, finally, convert into a Unix time.
    try:
        # Python 3.
        return dt.timestamp()
    except AttributeError:
        # Hack for supporting Python 2.
        return _datetime_to_timestamp(dt, tzinfo)


def _format_datetime_lazy(string):
    """Convenience string-processing function to adapt a string expressing a
    date correctly in our %Y-%m-%dT%H:%M:%S.%f format, adding the missing
    fields, if any.
    """
    if not 'T' in string:
        # This is a pure date---add all the missing fields.
        return '%sT00:00:00.0' % string
    if '.' in string:
        # We assume that this is a well-formed, complete datetime string.
        return string
    # If we are here it means that we have a partially formed datetime and we
    # need to fill in the blanks.
    fields = string.split('T')[1].split(':')
    assert len(fields) > 0
    return '%s%s.0' % (string, ':00' * (3 - len(fields)))


def string_to_unix(string, tzinfo=UTC_TZINFO, fmt=DATETIME_FMT):
    """Convert a string expressing time and date to a Unix time.

    Args
    ----
    string : string
        The input datetime string.

    tzinfo : a datetime.timezone instance or None
        The timezone info (use None for local time).

    fmt : string
        An optional format specifier for non standard input string

    Returns
    -------
    float
        The Unix time corresponding to the input string.
    """
    return _timestamp(string, tzinfo, fmt=fmt)


def string_to_unix_utc(string, fmt=DATETIME_FMT):
    """Convert a string expressing a UTC time and date to a Unix time.

    Args
    ----
    string : string
        The input datetime string.

    fmt : string
        An optional format specifier for non standard input string

    Returns
    -------
    float
        The Unix time corresponding to the input string.
    """
    return _timestamp(string, UTC_TZINFO, fmt=fmt)


def string_to_unix_local(string, fmt=DATETIME_FMT):
    """Convert a string expressing a local time and date to a Unix time.

    Args
    ----
    string : string
        The input datetime string.

    fmt : string
        An optional format specifier for non standard input string

    Returns
    -------
    float
        The Unix time corresponding to the input string.
    """
    return _timestamp(string, None, fmt=fmt)


def string_to_met(string, tzinfo=UTC_TZINFO, fmt=DATETIME_FMT):
    """Convert a string expressing time and date to a MET.

    Args
    ----
    string : string
        The input datetime string.

    tzinfo : a datetime.timezone instance or None
        The timezone info (use None for local time).

    fmt : string
        An optional format specifier for non standard input string

    Returns
    -------
    float
        The mission elapsed time corresponding to the input string.
    """
    return unix_to_met(string_to_unix(string, tzinfo, fmt=fmt))


def string_to_met_utc(string, lazy=False, fmt=DATETIME_FMT):
    """Convert a string expressing a UTC time and date to a MET.

    Args
    ----
    string : string
        The input datetime string.

    lazy : bool
        Flag to attempt and auto-fix partially formed datetime strings.

    fmt : string
        An optional format specifier for non standard input string

    Returns
    -------
    float
        The mission elapsed time corresponding to the input string.
    """
    if lazy:
        string = _format_datetime_lazy(string)
    return string_to_met(string, UTC_TZINFO, fmt=fmt)


LAUNCH_MET = string_to_met_utc(LAUNCH_DATETIME)


def string_to_met_local(string, fmt=DATETIME_FMT):
    """Convert a string expressing a local time and date to a MET.

    Args
    ----
    string : string
        The input datetime string.

    fmt : string
        An optional format specifier for non standard input string

    Returns
    -------
    float
        The mission elapsed time corresponding to the input string.
    """
    return string_to_met(string, None, fmt=fmt)


def current_time():
    """Return the current unix time.

    Returns
    -------
    float
        The current Unix time.
    """
    return time.time()


def current_met():
    """Return the current mission elapsed time.

    Returns
    -------
    float
        The current mission elapsed time.
    """
    return unix_to_met(current_time())


def current_datetime_string(tzinfo=UTC_TZINFO, fmt=DATETIME_FMT):
    """Return a string with the current date and time.

    Args
    ----
    tzinfo : a datetime.timezone instance or None
        The timezone info (use None for local time).

    Returns
    -------
    string
        The current date and time string.
    """
    return unix_to_string(current_time(), tzinfo, fmt)


def current_datetime_string_utc(fmt=DATETIME_FMT):
    """Return a string with the current UTC date and time.

    Args
    ----
    fmt : string
        An optional format specifier for non standard input string

    Returns
    -------
    string
        The current UTC date and time string.
    """
    return unix_to_string(current_time(), UTC_TZINFO, fmt)


def current_datetime_string_local(fmt=DATETIME_FMT):
    """Return a string with the current UTC date and time.

    Args
    ----
    fmt : string
        An optional format specifier for non standard input string

    Returns
    -------
    string
        The current local date and time string.
    """
    return unix_to_string(current_time(), None, fmt)


__SECONDS_IN_DAY = 86400.
__SECONDS_IN_YEAR = 31557600.


def seconds_to_days(seconds):
    """Convert seconds to days.
    """
    return seconds / __SECONDS_IN_DAY


def days_to_seconds(days):
    """Convert days to seconds.
    """
    return days * __SECONDS_IN_DAY


def seconds_to_years(seconds):
    """Convert seconds to days.
    """
    return seconds / __SECONDS_IN_YEAR


def years_to_seconds(years):
    """Convert years to seconds.
    """
    return years * __SECONDS_IN_YEAR


def met_to_mjd(met):
    """Convert a MET in the corresponding Modified Julian date.

    Args
    ----
    met : float
        The mission elapsed time.

    Returns
    -------
    float
        The Modified Julian Date.
    """
    return MISSION_START_MJD + seconds_to_days(met)


def mjd_to_met(mjd):
    """Convert a MJD in the corresponding Mission Elapsed Time.

    Args
    ----
    mjd : float
        The time ref to Modified Julian Day.

    Returns
    -------
    float
        The mission elapsed time.
    """
    met = mjd - MISSION_START_MJD - MISSION_START_MJDREFF
    return days_to_seconds(met)


MJD_OFFSET = 2400000.5


def met_to_jd(met):
    """Convert a MET in the corresponding Julian date.

    Args
    ----
    met : float
        The mission elapsed time.

    Returns
    -------
    float
        The Julian Date.
    """
    return met_to_mjd(met) + MJD_OFFSET


MET_TO_NUM_OFFSET = date2num(MISSION_START_DATETIME)


def met_to_num(met):
    """Convenience conversion factor to turn a MET into a number that matplotlib
    can interpret natively as a datetime.

    According to https://matplotlib.org/3.1.1/api/dates_api.html
    matplotlib represents dates using floating point numbers specifying the
    number of days since 0001-01-01 UTC, plus 1.
    """
    return MET_TO_NUM_OFFSET + seconds_to_days(met)

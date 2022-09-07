#!/usr/bin/env python
#
# Copyright (C) 2015--2018, the ixpeobssim team.
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

"""Collection of utilities for minimal code profiling.
"""

from functools import wraps
import time

try:
    import psutil
except ImportError:
    psutil = None

from ixpeobssim.utils.logging_ import logger


NO_PSUTIL_MSG = 'No system information available.'


def MB(size):
    """Convert from B to MB.
    """
    return size / 1048576.


def _pswrap(func):
    """
    """
    if psutil is not None:
        vm = psutil.virtual_memory()
        return func(vm)
    return -1


def psmem():
    """Return the virtual memory profile (in MB), as provided by psutil.
    """
    def func(vm):
        return MB(vm.total), MB(vm.used), MB(vm.free), MB(vm.available)
    return _pswrap(func)


def psfree():
    """Return the amount of free memory (in MB), as provided by psutil.
    """
    def func(vm):
        return MB(vm.free)
    return _pswrap(func)


def psavailable():
    """Return the amount of available memory (in MB), as provided by psutil.
    """
    def func(vm):
        return MB(vm.available)
    return _pswrap(func)


def psstatus():
    """Return a string representing the status of the memory.
    """
    if psutil is None:
        return NO_PSUTIL_MSG
    return 'Memory: %.2f MB, %.2f MB used, %.2f MB free %.2f MB available' %\
        psmem()



def timing(f):
    """Small decorator to time a generic function.
    """
    @wraps(f)
    def wrap(*args, **kwargs):
        start_time = time.time()
        result = f(*args, **kwargs)
        elapsed_time = time.time() - start_time
        logger.info('Running time for %s(): %.3f s', f.__name__, elapsed_time)
        return result
    return wrap


class xChrono:

    """Small chronometer class.

    A chronometer essentially measures the elapsed time since it has been
    started and is equipped to print itself to the standard output. (Note the
    chronometer is reset unpon the instantiation of a class object.)

    Examples
    --------
    >>> from ixpeobssim.utils.profile import xChrono
    >>> c = xChrono()
    >>> # ... do something.
    >>> print(c)
    """

    def __init__(self):
        """Constructor.
        """
        self.reset()

    def reset(self):
        """Reset the chronometer.
        """
        self.start_time = time.time()

    def __call__(self):
        """Return the elapsed time.
        """
        return time.time() - self.start_time

    def __str__(self):
        """ String formatting.
        """
        return '[t0 + %.3f s]' % self()



class xMemoryProfiler:

    """Small utility class to help profiling the allocated memory.
    """

    def __init__(self):
        """Constructor.
        """
        self.start_avail_mem = psavailable()
        self.last_avail_mem = self.start_avail_mem
        self.chrono = xChrono()
        if psutil is None:
            logger.info(NO_PSUTIL_MSG)

    @classmethod
    def available(self):
        """Return the available memory.
        """
        return psavailable()

    def __call__(self):
        """Implementation of class call.
        """
        self.last_avail_mem = self.available()
        return self.last_avail_mem, self.last_avail_mem - self.start_avail_mem

    def __str__(self):
        """Terminal formatting.
        """
        available, delta = self()
        return 'Available memory: %.2f MB, delta = %.2f MB %s' %\
            (available, delta, self.chrono)

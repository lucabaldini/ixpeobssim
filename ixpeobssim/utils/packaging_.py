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

import re

from ixpeobssim.utils.logging_ import abort


"""Verbatim from https://github.com/pypa/packaging/blob/master/packaging/version.py

This supposedly implements the indications in https://www.python.org/dev/peps/pep-0440/
"""
_VERSION_PATTERN = r"""
    v?
    (?:
        (?:(?P<epoch>[0-9]+)!)?                           # epoch
        (?P<release>[0-9]+(?:\.[0-9]+)*)                  # release segment
        (?P<pre>                                          # pre-release
            [-_\.]?
            (?P<pre_l>(a|b|c|rc|alpha|beta|pre|preview))
            [-_\.]?
            (?P<pre_n>[0-9]+)?
        )?
        (?P<post>                                         # post release
            (?:-(?P<post_n1>[0-9]+))
            |
            (?:
                [-_\.]?
                (?P<post_l>post|rev|r)
                [-_\.]?
                (?P<post_n2>[0-9]+)?
            )
        )?
        (?P<dev>                                          # dev release
            [-_\.]?
            (?P<dev_l>dev)
            [-_\.]?
            (?P<dev_n>[0-9]+)?
        )?
    )
    (?:\+(?P<local>[a-z0-9]+(?:[-_\.][a-z0-9]+)*))?       # local version
"""

_VERSION_REGEX = re.compile(r"^\s*" + _VERSION_PATTERN + r"\s*$", re.VERBOSE | re.IGNORECASE)



class xPackageVersion:

    """Small class encapsulating a package version.

    This was introduced when fixing issue #280, and the whole mechanism is
    inspired by the wonderful packaging Python package (since the use of the
    version information we are doing is fairly limited we didn't want to add
    yet another dependence on ixpeobssim).
    """

    def __init__(self, major, minor, patch=None):
        """Constructor.
        """
        self.major = major
        self.minor = minor
        self.patch = patch
        self._key = (major, minor, patch or 0)

    def __lt__(self, other):
        """Overloaded operator
        """
        return self._compare(other, lambda s, o: s < o)

    def __le__(self, other):
        """Overloaded operator
        """
        return self._compare(other, lambda s, o: s <= o)

    def __gt__(self, other):
        """Overloaded operator
        """
        return self._compare(other, lambda s, o: s > o)

    def __ge__(self, other):
        """Overloaded operator
        """
        return self._compare(other, lambda s, o: s >= o)

    def __eq__(self, other):
        """Overloaded operator
        """
        return self._compare(other, lambda s, o: s == o)

    def __ne__(self, other):
        """Overloaded operator
        """
        return self._compare(other, lambda s, o: s != o)

    def _compare(self, other, method):
        """
        """
        if isinstance(other, str):
            other = parse_version_string(other)
        return method(self._key, other._key)

    def __str__(self):
        """String formatting.
        """
        return 'Version: major=%d, minor=%d, patch=%s' %\
            (self.major, self.minor, self.patch)



def parse_version_string(version_string):
    """Small utility function to parse the version string of a generic package.

    Note that we only parse the major, minor and patch fields (i.e., the release
    segment) of the version string.
    """
    match = _VERSION_REGEX.search(version_string)
    if not match:
        abort('Cannot parse version string %s' % version_string)
    release = tuple(int(i) for i in match.group('release').split('.'))
    return xPackageVersion(*release)


def retrieve_version(package):
    """Retrieve the version for a generic package.
    """
    return parse_version_string(package.__version__)

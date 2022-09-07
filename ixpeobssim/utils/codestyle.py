#!/usr/bin/env python
#
# Copyright (C) 2015, the ixpeobssim team.                              *
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


"""Small fake module to illustrate the basic coding guidelines and the use
of docstrings for sphinx.

See `PEP 0008 <https://www.python.org/dev/peps/pep-0008/>`_ for more details.
(And note that the module name is short, descriptive, and all lowercase).
Also note we use the Napoleon sphinx extension with the Numpy style.
"""

# One import per line, all at the top of the module.
import math


"""The electron mass---this is a global constant, and all uppercase.
"""
ELECTRON_MASS = 0.510998910 # MeV


def square(x):
    """Return the square of x.

    Args
    ----
    x : number or array
        The input x value (can be either a numeric type or an array).

    Returns
    -------
    float
        The square of the input value

    Examples
    --------
    >>> x = 3.
    >>> y = square(x)
    >>> print(y)

    Note
    ----
    Yes, this is a pretty dumb function, but watch out for the use of doctrings.
    """
    return math.pow(x, 2.)


class xSampleClass:

    """An example class, doing nothing useful.

    Args
    ----
    name : string
        The instance name.
    description : string
        An optional description.

    Examples
    --------
    >>> from codestyle import xSampleClass
    >>> c1 = xSampleClass('hello')
    >>> c2 = xSampleClass('world!')
    >>> print(c1 + c2)
    """

    def __init__(self, name, description=None):
        """Constructor.

        Note there are no spaces around the `=` sign for the default parameter
        value.
        """
        self.name = name
        self.description = description

    def __add__(self, other):
        """Implementation of the `+` operator.

        This is just concatenating the string representation of the two
        objects being added.
        """
        return self.__class__('%s %s' % (self, other))

    def __str__(self):
        """String formatting.
        """
        return self.name


if __name__ == '__main__':
    c1 = xSampleClass('hello')
    c2 = xSampleClass('world!')
    print(c1 + c2)

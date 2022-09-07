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

"""Geometry module.
"""

from __future__ import print_function, division

import numpy


class xPoint:

    """Class representing a point in the three-dimensional space.
    """

    def __init__(self, x, y, z):
        """Constructor.
        """
        self.x = x
        self.y = y
        self.z = z

    def norm(self):
        """Return the norm of the three-dimensional vector corresponding to the
        point.
        """
        return numpy.sqrt(self.x**2. + self.y**2. + self.z**2.)

    def distance_to(self, other):
        """Return the distance to another point.
        """
        return (self - other).norm()

    def move(self, length, ray):
        """Move the point by a given length along a given ray.
        """
        x = self.x + length * ray.xdir
        y = self.y + length * ray.ydir
        z = self.z + length * ray.zdir
        return self.__class__(x, y, z)

    @classmethod
    def unphysical_point(cls):
        """Return an unphysical point.
        """
        return cls(numpy.nan, numpy.nan, numpy.nan)

    def unphysical(self):
        """Return True if the point is unphysical.
        """
        return numpy.nan in [self.x, self.y, self.z]

    def __add__(self, other):
        """Operator overload.
        """
        return self.__class__(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        """Operator overload.
        """
        return self.__class__(self.x - other.x, self.y - other.y, self.z - other.z)

    def __eq__(self, other):
        """Operator overload.
        """
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __str__(self):
        """String formatting.
        """
        args = self.__class__.__name__, self.x, self.y, self.z
        return '%s(%.6f, %.6f, %.6f)' % args



class xLine:

    """Class representing a line.
    """

    def __init__(self, begin, end):
        """Constructor.
        """
        self.begin = begin
        self.end = end

    def length(self):
        """Return the length of the line.
        """
        return (self.end - self.begin).norm()

    def __str__(self):
        """String formatting.
        """
        return '%s--%s' % (self.begin, self.end)



class xRay:

    """Class representing a ray in the three-dimensional space.
    """

    def __init__(self, origin, theta, phi):
        """Constructor.
        """
        self.origin = origin
        self.theta = theta
        self.phi = phi
        # Calculate and cache the cosine directors for later use.
        st = numpy.sin(numpy.radians(self.theta))
        ct = numpy.cos(numpy.radians(self.theta))
        sp = numpy.sin(numpy.radians(self.phi))
        cp = numpy.cos(numpy.radians(self.phi))
        self.xdir = st * cp
        self.ydir = st * sp
        self.zdir = ct

    def xintersect(self, x):
        """Return the point where the ray intersects the plane at a given x.
        """
        if x == self.origin.x:
            return self.origin
        if self.xdir == 0.:
            return xPoint.unphysical_point()
        dx = (x - self.origin.x) / self.xdir
        return xPoint(x, self.origin.y + dx * self.ydir, self.origin.z + dx * self.zdir)

    def yintersect(self, y):
        """Return the point where the ray intersects the plane at a given y.
        """
        if y == self.origin.y:
            return self.origin
        if self.ydir == 0.:
            return xPoint.unphysical_point()
        dy = (y - self.origin.y) / self.ydir
        return xPoint(self.origin.x + dy * self.xdir, y, self.origin.z + dy * self.zdir)

    def zintersect(self, z):
        """Return the point where the ray intersects the plane at a given z.
        """
        if z == self.origin.z:
            return self.origin
        if self.zdir == 0.:
            return xPoint.unphysical_point()
        dz = (z - self.origin.z) / self.zdir
        return xPoint(self.origin.x + dz * self.xdir, self.origin.y + dz * self.ydir, z)

    def __str__(self):
        """String formatting.
        """
        args = self.__class__.__name__, self.origin, self.theta, self.phi
        return '%s: %s -> (%.3f deg, %.3f deg)' % args

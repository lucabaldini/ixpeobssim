#!/usr/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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

import numpy

from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.evt.animate import xMovingCircle
from ixpeobssim.utils.matplotlib_ import plt


"""Configuration file for a Cassiopea A animation.

This is describing a moving circular region of interest that can be used in
conjunction with a Cas A map to create an audio-visual animation.

.. warning::

   The code should be cleaned up and streamlined, but hopefully you get the
   main points :-)
"""

# Basic settings.
kxy = 2
kr = 1
t1 = 50.
t2 = 60.
t3 = 62.
t4 = 65.

# Definition of the initial part of the spatial path.
xy = (2., -2.), (-1., -1.), (-2.5, 0.), (-1.5, 1.5), (1.75, 2.25), (4., 1.), (3., -0.5), (0.5, 0.)
x = numpy.array([pos[0] for pos in xy])
y = numpy.array([pos[1] for pos in xy])

# Calculate the appropriate times for the roi to move at a constant velocity
# along the path.
oversample = 100
t = numpy.linspace(0., t1, len(xy))
sx = xInterpolatedUnivariateSpline(t, x, k=kxy)
sy = xInterpolatedUnivariateSpline(t, y, k=kxy)
t = numpy.linspace(0., t1, oversample * len(xy))
gridx, gridy = sx(t), sy(t)
dx, dy = numpy.diff(gridx), numpy.diff(gridy)
dr = numpy.sqrt(dx**2. + dy**2.)
r = numpy.array([0.] + [dr[0:(i + 1) * oversample].sum() for i in range(len(xy) - 1)])
t = r * t1 / r[-1]

#plt.plot(gridx, gridy)
#plt.plot(x, y, 'o')
#plt.gca().set_aspect('equal')
#plt.show()

roi_positions = list(zip(t, x, y)) + [(t2, 0., 0.), (t3, 0., 0.), (t4, 0., 0.)]
roi_radius = [25.] * len(xy) + [250., 0.1, 0.1]
ROI = xMovingCircle(roi_positions, roi_radius, kxy=2)

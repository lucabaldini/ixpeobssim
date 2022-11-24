#!/urs/bin/env python
#
# Copyright (C) 2018--2022, the ixpeobssim team.
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

"""Basic GPD-related constants.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.instrument.du import du_rotation_angle
from ixpeobssim.utils.math_ import modulo_2pi

# pylint: disable=invalid-name

def rectangle_area(half_side_x, half_side_y):
    """Small convenience function to calulate the area of a rectangle, given the
    half sides along the two coordinates.

    Arguments
    ---------
    half_side_x : float
        The half side of the rectangle along the x coordinate (in mm).

    half_side_y : float
        The half side of the rectangle along the y coordinate (in mm).
    """
    return 4. * half_side_x * half_side_y


# Physical side of the readout chip on the two orthogonal directions...
GPD_PHYSICAL_HALF_SIDE_X = 7.4875
GPD_PHYSICAL_HALF_SIDE_Y = 7.599
# ... and fellow derived quantities.
GPD_PHYSICAL_AREA = rectangle_area(GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_Y)
GPD_PHYSICAL_MAX_RADIUS = numpy.sqrt(GPD_PHYSICAL_HALF_SIDE_X**2. +\
    GPD_PHYSICAL_HALF_SIDE_Y**2.)

# And not the default values for the fiducial half side in detector coordinates.
# Note the fiducial rectangle has changed along the way due to a small bug in the
# processing code and, unlike the physical dimensions of the readout chip, is not
# guaranteed to be the same for all observations.
GPD_DEFAULT_FIDUCIAL_HALF_SIDE_X = 6.600
GPD_DEFAULT_FIDUCIAL_HALF_SIDE_Y = 6.800

# GEM manufacturing parametes
LASER_ETCHING_PITCH = 1.800
NUM_LASER_SWEEPS = 8


def fiducial_area(half_side_x=GPD_DEFAULT_FIDUCIAL_HALF_SIDE_X,
    half_side_y=GPD_DEFAULT_FIDUCIAL_HALF_SIDE_Y):
    """Return the area of the fiducial rectangle.

    This is essentially calling the base rectangle_area() function, with the
    proper default values for the fiducial cut.

    Arguments
    ---------
    half_side_x : float
        The half side of the fiducial rectangle along the x coordinate in mm.

    half_side_y : float
        The half side of the fiducial rectangle along the y coordinate in mm.
    """
    return rectangle_area(half_side_x, half_side_y)


def gpd_map_binning(half_side_x, half_side_y, num_bins_x, num_bins_y=None):
    """Return a numpy array with an appropriate binning for a bidimensional
    map over the GPD area.

    Note this function was refactored to adapt to a generic rectangle in response
    to https://github.com/lucabaldini/ixpeobssim/issues/668, so this now
    returns two independent arrays representing the binng on the x and y axes.

    Arguments
    ---------
    half_side_x : float
        The half side of the histogram binning along the x coordinate (in mm).

    half_side_y : float
        The half side of the histogram binning along the y coordinate (in mm).
    """
    if num_bins_y is None:
        num_bins_y = num_bins_x
    return numpy.linspace(-half_side_x, half_side_x, num_bins_x + 1),\
        numpy.linspace(-half_side_y, half_side_y, num_bins_y + 1)


def within_fiducial_rectangle(x, y, half_side_x=GPD_DEFAULT_FIDUCIAL_HALF_SIDE_X,
    half_side_y=GPD_DEFAULT_FIDUCIAL_HALF_SIDE_Y):
    """Return wheter an (x, y) position in mm is within a fiducial rectangle of
    given half-sides on the two coordinates.

    Arguments
    ---------
    x : float or array
        The x position or array of x positions in mm

    y : float or array
        The y position or array of y positions in mm

    half_side_x : float
        The half side of the fiducial rectangle along the x coordinate (in mm).

    half_side_y : float
        The half side of the fiducial rectangle along the y coordinate (in mm).
    """
    return numpy.logical_and(abs(x) <= half_side_x, abs(y) <= half_side_y)


def within_gpd_physical_area(x, y):
    """Return wheter an (x, y) position in mm is within the GPD physical area.

    Arguments
    ---------
    x : float or array
        The x position or array of x positions in mm

    y : float or array
        The y position or array of y positions in mm
    """
    return within_fiducial_rectangle(x, y, GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_Y)


def rotate_detxy(x, y, du_id, roll_angle=0., inverse=False):
    """Convert the (x, y) position from the focal plane reference frame to the
    the gpd reference frame based on the DU id and roll angle.

    .. warning::
       This would be more consistent with the rest of the code base if the direct
       and inverse rotations were implemented in two different functions rather
       than one function with an additional argument.

    Arguments
    ---------
    x : float or array
        The x position or array of x positions in mm

    y : float or array
        The y position or array of y positions in mm

    du_id : int
        The detector unit hosting the GPD

    roll_angle : float
        The telescope roll angle in decimal degrees
    """
    rotation_angle = du_rotation_angle(du_id, roll_angle)
    if inverse:
        rotation_angle = -rotation_angle
    sinr = numpy.sin(rotation_angle)
    cosr = numpy.cos(rotation_angle)
    rx = cosr * x - sinr * y
    ry = sinr * x + cosr * y
    return rx, ry


def phi_to_detphi(phi, du_id, roll_angle):
    """Convert the azimuthal angle from the focal plane (i.e., sky) reference
    frame to the GPD reference.

    Arguments
    ---------
    phi : float or array
        The azimuthal angle in the sky reference frame in radians

    du_id : int
        The detector unit hosting the GPD

    roll_angle : float
        The telescope roll angle in decimal degrees
    """
    rotation_angle = du_rotation_angle(du_id, roll_angle)
    detphi = phi + rotation_angle
    return modulo_2pi(detphi)


def detphi_to_phi(detphi, du_id, roll_angle):
    """Convert the phi angle from the GPD reference frame to the focal plane
    (i.e., sky) reference.

    Arguments
    ---------
    detphi : float or array
        The azimuthal angle in the GPD reference frale in radians

    du_id : int
        The detector unit hosting the GPD

    roll_angle : float
        The telescope roll angle in decimal degrees
    """
    rotation_angle = -du_rotation_angle(du_id, roll_angle)
    phi = detphi + rotation_angle
    return modulo_2pi(phi)

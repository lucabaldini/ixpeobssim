#!/urs/bin/env python
#
# Copyright (C) 2018--2019, the ixpeobssim team.
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

from ixpeobssim.instrument.du import du_rotation_angle
from ixpeobssim.utils.math_ import modulo_2pi


PHYSICAL_HALF_SIZE = 7.500
FIDUCIAL_HALF_SIZE = 7.350
LASER_ETCHING_PITCH = 1.800
NUM_LASER_SWEEPS = 8
FIDUCIAL_AREA = (2. * FIDUCIAL_HALF_SIZE)**2.


def gpd_map_binning(num_bins):
    """Return a numpy array with an appropriate binning for a bidimensional
    map over the GPD area.
    """
    return numpy.linspace(-FIDUCIAL_HALF_SIZE, FIDUCIAL_HALF_SIZE, num_bins + 1)


def within_fiducial_area(x, y):
    """Return wheter an (x, y) position in mm is within the fiducial active
    area of the GPD readout ASIC.

    Arguments
    ---------
    x : float or array
        The x position or array of x positions in mm

    y : float or array
        The y position or array of y positions in mm
    """
    return (abs(x) <= FIDUCIAL_HALF_SIZE) * (abs(y) <= FIDUCIAL_HALF_SIZE)


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

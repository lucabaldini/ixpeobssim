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

from __future__ import print_function, division

import numpy

from ixpeobssim.instrument import DU_IDS
from ixpeobssim.utils.math_ import modulo_2pi


# Rotation angle for the three detector units. These angles are measured in the
# focal plane reference frame and correspond to the direction of the x-axes
# (phi = 0) of the three gpd reference frames.
# Derived from drawing number 2506942 (Rev. C).
__DU_ROTATION_ANGLE = {
    DU_IDS[0]: numpy.radians(109.), # 90 + 19
    DU_IDS[1]: numpy.radians(229.), # 120 + 90 + 19
    DU_IDS[2]: numpy.radians(349.)  # 240 + 90 + 19
}


# Mapping of the detector units from logical to physical space.
__DU_MODEL = {1: 'DU_FM2', 2: 'DU_FM3', 3: 'DU_FM4'}


def du_rotation_angle(du_id, roll_angle=0., modulo=True):
    """Return the rotation angle for a given DU (modulo 2 pi).
    """
    if not du_id in DU_IDS:
        raise RuntimeError('Invalid DU ID (%d)' % du_id)
    angle = __DU_ROTATION_ANGLE[du_id] + numpy.radians(roll_angle)
    if modulo:
        angle = modulo_2pi(angle)
    return angle


def du_logical_name(du_id):
    """Return the logical name for a given DU (i.e., the content of the DETNAME
    header keyword).
    """
    if not du_id in DU_IDS:
        raise RuntimeError('Invalid DU ID (%d)' % du_id)
    return 'DU%d' % du_id


def du_physical_name(du_id):
    """Return the physical name for a given DU (i.e., the content of the DET_ID
    header keyword).
    """
    if not du_id in DU_IDS:
        raise RuntimeError('Invalid DU ID (%d)' % du_id)
    return __DU_MODEL[du_id]


def det_name_to_du_id(det_name):
    """Convert the logical name of a given DU to the actual DU ID.
    """
    if not det_name.startswith('DU'):
        raise RuntimeError('Invalid detector name ("%s")' % det_name)
    return int(det_name[-1])

#!/usr/bin/env python
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


"""Some basic detector things.
"""
DU_IDS = [1, 2, 3]

NUM_DETECTOR_UNITS = len(DU_IDS)

COMBINED_DU_ID = 123


def du_suffix(du_id):
    """Return a suffix to identify a specific DU in a file name.
    """
    return 'du%d' % du_id


def du_id_is_valid(du_id):
    """Return True if the argument is a valid DU ID.
    """
    return du_id in DU_IDS + [COMBINED_DU_ID]

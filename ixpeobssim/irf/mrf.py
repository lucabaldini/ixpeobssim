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

"""Modulation respose function.
"""

from __future__ import print_function, division


from ixpeobssim.irf.base import xSpecRespBase


# pylint: disable=invalid-name


class xModulationResponse(xSpecRespBase):

    """Class describing the modulation response, i.e., the product of the
    effective area times the modulation factor.
    """

    Y_UNITS = 'cm$^2$'
    Y_LABEL = 'Modulation response function [%s]' % Y_UNITS

    def __init__(self, file_path):
        """Constructor.
        """
        xSpecRespBase.__init__(self, file_path, 'mrf')

    def plot(self):
        """Plot the modulation response.
        """
        # pylint: disable=arguments-differ
        self.plot_base(logy=True)

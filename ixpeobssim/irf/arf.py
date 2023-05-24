# Copyright (C) 2015--2022, the ixpeobssim team.
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

"""Read interface to effective area files.
"""

from __future__ import print_function, division

from ixpeobssim.irf.base import xSpecRespBase


# pylint: disable=invalid-name


class xEffectiveArea(xSpecRespBase):

    WEIGHTING_SCHEME_HEADER_KEYWORD = 'STOKESWT'
    VALID_WEIGHTING_SCHEMES = (None, 'UNWEIGHTED', 'NEFF', 'SIMPLE')

    """Class describing the on-axis effective area.

    Arguments
    ---------
    file_path : str
        The path to the .arf FITS file containing the effective area table.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xSpecRespBase.__init__(self, file_path, 'arf')

    def weighting_scheme(self):
        """Return the weighting scheme used to calculate the effective area.

        See https://bitbucket.org/ixpesw/ixpeobssim/issues/573/
        for more information about why we do need to keep track of this.

        The basic purpose that this hook serves is to be able to make sure that
        weights are used consistently across binned data products and
        response files, e.g., we don't want for people to be able to create
        weighted polarization cubes with the effective area correction unless
        the arf file is created with the 'SIMPLE' prescription.

        For backward compatibility this is returning None if the response file
        is missing the `STOKESWT` keyword in the `SPECRESP` extension.
        If this is the case, the user should not make assumptions about the
        actual weighting scheme used.
        """
        header = self.hdu_list['SPECRESP'].header
        scheme = header.get(self.WEIGHTING_SCHEME_HEADER_KEYWORD, None)
        if not scheme in self.VALID_WEIGHTING_SCHEMES:
            raise RuntimeError('Invalid arf weighting scheme: %s' % scheme)
        return scheme

    def plot(self):
        """Plot the effective area.
        """
        # pylint: disable=arguments-differ
        self.plot_base(logy=True)



class xTowEffectiveArea(xSpecRespBase):

    """Class describing the on-axis effective area at the top of the window.

    Arguments
    ---------
    file_path : str
        The path to the FITS file containing the effective area table.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xSpecRespBase.__init__(self, file_path, 'fits')

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

"""Base classes for all the read interfaces to response files.

The classes in this module represent the base of the hyerarchy for all the
classes descibing response files, and particularly:

* xResponseBase acts as a base class for *all* the response files;
* xSpecRespBase acts as a base class for the (on-axis) effective area, the
  modulation factor and the modulation respose function.

The latter (i.e., the classes inheriting from xSpecRespBase) all have in common
being one-dimensional functions of the (true) energy, and have a SPECRESP
extension containing the response in bins of energy.

The xSpecRespBase class is also equipped to allow for systematic errors on the
response values, according to the proper OGIP standard.
"""

from __future__ import print_function, division

import os

import numpy

from ixpeobssim.core.fitsio import read_hdu_list_in_memory
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.instrument.du import det_name_to_du_id
from ixpeobssim.irf.ebounds import ENERGY_MIN, ENERGY_MAX
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, last_line_color, setup_gca, \
    residual_plot, labeled_marker


# pylint: disable=invalid-name, no-member


class xResponseBase:

    """Base class for reading response data from FITS files.

    Arguments
    ---------
    file_path : str
        The path to the input FITS file.

    extension : str
        The extension of the input FITS file (typically fits, arf or rmf)
    """

    def __init__(self, file_path, extension):
        """Constuctor.
        """
        self.file_path = file_path
        self.base_name = os.path.basename(file_path)
        self.hdu_list = read_hdu_list_in_memory(file_path, extension)
        self.primary_header = self.hdu_list['PRIMARY'].header
        self.irf_type = self.primary_header['IRFTYPE']
        try:
            self.du_id = det_name_to_du_id(self.primary_header['DETNAM'])
        except AssertionError:
            # Hack to handle the case where, for old response files, we shipped
            # the combined response for the three detector units.
            self.du_id = None

    def header_comments(self):
        """Return the content of the COMMENT keyword of the primary header, split
        into lines.

        This is essentially only used in ixpeobssim.evt.ixpesim to reverse-engineer
        the value of the GPD pressure for the purpose of building the telescope
        response at the top of the Be window.
        """
        try:
            return str(self.primary_header['COMMENT']).split('\n')
        except KeyError:
            return None

    def field(self, col_name, ext_index=1):
        """Return a view over a column of the specified extension as an array.

        Since most of the response files have exactly one relevant (binary table)
        extension, having the extension index defaulting to one is a handy
        shortcut to retrieve the column data in a general fashion.

        Args
        ----
        col_name : str
            The name of the target column in the binary table extension.

        ext_index : int
            The index of the target binary table extension.
        """
        return self.hdu_list[ext_index].data.field(col_name)

    def __str__(self):
        """String representation.
        """
        return 'Response file %s %s' % (self.base_name, self.primary_header)



class xSpecRespBase(xResponseBase, xInterpolatedUnivariateSpline):

    """Derived class describing a spectral response, i.e., effective area,
    modulation response function, or modulation factor.
    """

    Y_UNITS = 'cm$^2$'
    Y_LABEL = 'On-axis effective area [%s]' % Y_UNITS

    def __init__(self, file_path, extension, k=2):
        """Overloaded constructor
        """
        xResponseBase.__init__(self, file_path, extension)
        energy = 0.5 * (self.field('ENERG_LO') + self.field('ENERG_HI'))
        try:
            resp = self.field('SPECRESP')
        except KeyError:
            # Horrible hack to support pre-12 modulation factors, see
            # https://bitbucket.org/ixpesw/ixpeobssim/issues/339
            resp = self.field('MODFRESP')
        argmax = numpy.argmax(resp)
        self.peak_energy = energy[argmax]
        self.peak_value = resp[argmax]
        fmt = dict(xlabel='Energy [keV]', ylabel=self.Y_LABEL, k=k)
        xInterpolatedUnivariateSpline.__init__(self, energy, resp, **fmt)
        try:
            fmt = dict(xlabel='Energy [keV]', ylabel='Relative systematic error', k=k)
            self.sys_min = xInterpolatedUnivariateSpline(energy, self.field('SYS_MIN'), **fmt)
            self.sys_max = xInterpolatedUnivariateSpline(energy, self.field('SYS_MAX'), **fmt)
        except KeyError:
            self.sys_min = None
            self.sys_max = None

    def has_sys_errors(self):
        """Return True if the object instance has both negative and positive
        systematic errors defined.
        """
        return self.sys_min is not None and self.sys_max is not None

    def sys_envelope(self):
        """Return the envelope of the systematic uncertainties.
        """
        return self.y * (1. - self.sys_min(self.x)), self.y * (1. + self.sys_max(self.x))

    def plot_sys_envelope(self):
        """Plot the envelope of the systematic errors.
        """
        if not self.has_sys_errors():
            logger.warning('Cannot plot systematic bounds...')
            return
        _min, _max = self.sys_envelope()
        fmt = dict(linestyle='dashed', color=last_line_color())
        plt.plot(self.x, _min, label='$\\pm 1 \\sigma$ confidence interval', **fmt)
        plt.plot(self.x, _max, **fmt)
        plt.legend()

    def plot_sys_errors(self):
        """Plot the fractional systematic errors.
        """
        if not self.has_sys_errors():
            logger.warning('Cannot plot systematic bounds...')
        plt.plot(self.x, numpy.zeros(self.x.shape))
        fmt = dict(linestyle='dashed', color=last_line_color())
        plt.plot(self.x, -self.sys_min(self.x), **fmt)
        plt.plot(self.x, self.sys_max(self.x), **fmt)
        setup_gca(xmin=ENERGY_MIN, xmax=ENERGY_MAX, grids=True, xlabel='Energy [keV]',
            ylabel='Sys. error')

    def plot_base(self, **kwargs):
        """Plot the spectral response.
        """
        if self.has_sys_errors():
            _, ax2 = residual_plot(self.base_name)
        else:
            plt.figure(self.base_name)
        xInterpolatedUnivariateSpline.plot(self, label=self.base_name, **kwargs)
        label = '%.2f %s @ %.2f keV' % (self.peak_value, self.Y_UNITS, self.peak_energy)
        labeled_marker(self.peak_energy, self.peak_value, label, dy=2.5)
        setup_gca(xmin=ENERGY_MIN, xmax=ENERGY_MAX, grids=True, legend=True)
        if self.has_sys_errors():
            self.plot_sys_envelope()
            plt.sca(ax2)
            self.plot_sys_errors()

#!/usr/bin/env python
#
# Copyright (C) 2017--2019, the ixpeobssim team.
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


"""Small convenience module for interfacing with the XCOM database.

https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html
"""

import numpy
import os

from ixpeobssim import IXPEOBSSIM_IRFGEN
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.utils.matplotlib_ import plt


IXPEOBSSIM_XCOM_DATA = os.path.join(IXPEOBSSIM_IRFGEN, 'data', 'xcom')

__CACHE = {}


def create_energy_grid(emin=0.001, emax=0.015, num_points=100):
    """Create a text file with a list of additional energies to query the
    XCOM databse, to be passed to the XCOM web interface.

    By default this is writing out a grid of logarithmically-spaced energy
    values between 1 and 15 keV---one energy per line, in MeV.
    """
    file_path = os.path.join(IXPEOBSSIM_XCOM_DATA, 'egrid.txt')
    logger.info('Writing XCOM energy grid to %s...' % file_path)
    with open(file_path, 'w') as output_file:
        for energy in numpy.logspace(numpy.log10(emin), numpy.log10(emax),
                                     num_points):
            output_file.write('%.6f\n' % energy)
    logger.info('Done.')


class xCrossSectionTable:

    def __init__(self, identifier):
        """Constructor.

        Note we are converting Mev to keV in place.
        """
        self.identifier = identifier
        file_name = '%s.txt' % identifier.lower()
        file_path = os.path.join(IXPEOBSSIM_XCOM_DATA, file_name)
        logger.info('Parsing XCOM data file %s...' % file_path)
        self.energy, self.coherent, self.incoherent, self.photoelectric = \
            numpy.loadtxt(file_path, unpack=True)
        self.energy *= 1000.
        self.total = self.coherent + self.incoherent + self.photoelectric
        self.photoelectric_spline = self.__spline(self.photoelectric)

    def __spline(self, values):
        """Return a spline with a given cross section as a function of the
        energy.
        """
        fmt = dict(xlabel='Energy [keV]',
                   ylabel='Cross section [cm$^2$ g$^{-1}$]')
        return xInterpolatedUnivariateSplineLinear(self.energy, values, **fmt)

    def transparency(self, energy, thickness, density):
        """Return the transparency of a slab of material with a given density,
        evaluated on a grid of energy values.

        Mind this is based on the phototelectric cross section only
        (i.e., coherent and incoherent scattering are not included.)

        Args
        ----
        energy : array-like
            The array of energy values to evaluate the efficiency.

        thickness : float
            The thickness of the slab in cm.

        density : float
            The density of the material in g cm^{-3}
        """
        return numpy.exp(-density * self.photoelectric_spline(energy) * thickness)

    def photoabsorption_efficiency(self, energy, thickness, density):
        """Return the photoabsorption efficiency for a slab of material of
        a given thickness and density, evaluated on a grid of energy values.

        Args
        ----
        energy : array-like
            The array of energy values to evaluate the efficiency.

        thickness : float
            The thickness of the slab in cm.

        density : float
            The density of the material in g cm^{-3}
        """
        return 1. - self.transparency(energy, thickness, density)

    def plot(self):
        """Plot the cross section.
        """
        plt.figure(self.identifier)
        fmt = dict(logx=True, logy=True)
        self.__spline(self.coherent).plot(label='Coherent', **fmt)
        self.__spline(self.incoherent).plot(label='Incoherent', **fmt)
        self.__spline(self.photoelectric).plot(label='Photoelectric', **fmt)
        self.__spline(self.total).plot(label='Total', **fmt)
        plt.axis([self.energy.min(), self.energy.max(), None, None])
        plt.legend()

    def __str__(self):
        """Terminal formatting.
        """
        data = [self.energy, self.coherent, self.incoherent,\
                self.photoelectric, self.total]
        header = 'Energy [keV]   Coherent     Incoherent   Photoel.     Total'
        return '%s\n%s' % (header, numpy.around(numpy.vstack(data).T, 3))



def load_xsection_data(identifier):
    """Load the cross-section data for a given element or compound.
    """
    if identifier not in __CACHE:
        __CACHE[identifier] = xCrossSectionTable(identifier)
    return __CACHE[identifier]

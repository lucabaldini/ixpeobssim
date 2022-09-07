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

from __future__ import print_function, division


"""Small convenience module for interfacing with the astar database.

https://physics.nist.gov/PhysRefData/Star/Text/ASTAR.html
"""

import numpy
import os

from ixpeobssim import IXPEOBSSIM_IRFGEN
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.irfgen.constants import H_MASS, C_MASS, O_MASS, DME_MASS
from ixpeobssim.irfgen.constants import dme_density


IXPEOBSSIM_ASTAR_DATA = os.path.join(IXPEOBSSIM_IRFGEN, 'data', 'astar')


class xAlphaStoppingPowerTable:

    """Basic interface to the output files of the ASTAR database.
    """

    def __init__(self, identifier, density=None):
        """Constructor.

        There has been a lot of back and forth on this class, essentially
        because we want to code things in a sensible way for elements and at the
        same time we want to make it practical to handle compounds (e.g., DME).

        A few noticeable things:

        * here we are using MeV and cm (instead of keV and mm) throughout, as
        these are the "natural" units in this business;
        * we support an optional density argument in the constructor, that
        allows to handle solids in a natural fashion; this implies some
        bookkeeping with the units, as the numbers have a different meaning
        depending of whether they are normalized to the density or not.
        """
        self.identifier = identifier
        self.density = density
        # Read the underlying data file.
        file_name = '%s.txt' % identifier.lower()
        file_path = os.path.join(IXPEOBSSIM_ASTAR_DATA, file_name)
        logger.info('Parsing ASTAR data file %s...' % file_path)
        self._energy, _, _, self._stopping_power, self._csda_range, _, _ = \
            numpy.loadtxt(file_path, unpack=True)
        self.process_data()

    def process_data(self):
        """Do all the post-processing of the underlying data.
        """
        # Some book-keeping.
        energy_label = '$\\alpha$ kinetic energy [MeV]'
        stopping_power_label = 'Stopping power [%s]'
        range_label = 'CSDA range [%s]'
        if self.density is None:
            stopping_power_label = stopping_power_label % 'MeV cm$^2$ g$^{-1}$'
            range_label = range_label % 'g cm$^{-2}$'
        else:
            stopping_power_label = stopping_power_label % 'MeV cm$^{-1}$'
            range_label = range_label % 'cm'
            # If we define the density of the material, we need to rescale the
            # relevant quantities.
            self._stopping_power *= self.density
            self._csda_range /= self.density
        # Build a spline with the stopping power as a function of the energy.
        args = self._energy, self._stopping_power
        fmt = dict(xlabel=energy_label, ylabel=stopping_power_label)
        self.stopping_power = xInterpolatedUnivariateSpline(*args, **fmt)
        # Calculate the inverse stopping power spline, which will be handy for
        # later use.
        args = self._energy, 1. / self._stopping_power
        fmt = dict(xlabel=energy_label, ylabel='Inverse stopping power')
        self.inverse_stopping_power = xInterpolatedUnivariateSpline(*args, **fmt)
        # Now to the CSDA range spline.
        args = self._energy, self._csda_range
        fmt = dict(xlabel=energy_label, ylabel=range_label)
        self.csda_range = xInterpolatedUnivariateSpline(*args, **fmt)

    def bragg_curve(self, energy=10., energy_step=0.01):
        """Return a spline representing the Bragg curve for the element or
        compound.
        """
        depth = 0.
        x = []
        y = []
        while energy > 0.:
            stopping_power = self.stopping_power(energy)
            x.append(depth)
            if energy < energy_step:
                y.append(0.)
                break
            y.append(stopping_power)
            energy -= energy_step
            depth += energy_step / stopping_power
        x = numpy.array(x)
        y = numpy.array(y)
        fmt = dict(xlabel='Path Length [cm]', ylabel=self.stopping_power.ylabel)
        return xInterpolatedUnivariateSpline(x, y, **fmt)

    def _energy_profile(self, energy=10., energy_step=0.01):
        """Return a spline with the longitudinal energy profile of the alpha at
        a given energy as a function of the depth in the material.

        Note that for this to make sense the density has to be defined.
        """
        assert self.density is not None
        depth = 0.
        range_ = self.csda_range(energy)
        x = [depth]
        y = [energy]
        while energy > 0.:
            stopping_power = self.stopping_power(energy)
            energy -= energy_step
            depth += energy_step / stopping_power
            if energy < energy_step:
                # Add the last step.
                x.append(range_)
                y.append(0.)
                break
            x.append(depth)
            y.append(energy)
        x = numpy.array(x)
        y = numpy.array(y)
        fmt = dict(xlabel='Path Length [cm]', ylabel=self.stopping_power.xlabel)
        # Note that we set the ext spline parameter to 'const' so that the
        # energy profile returns 0. beyond the particle range.
        return xInterpolatedUnivariateSpline(x, y, ext='const', **fmt)

    def energy_loss(self, max_energy=10., energy_step=0.01):
        """
        """
        # Build the energy profile.
        _profile = self._energy_profile(max_energy, energy_step)
        # Invert the energy profile (mind you have to flip the arrays, as
        # y is decreasing).
        # Note we pass the axis=0 argument explicitly in order for this to
        # work with numpy versions prior to 1.15
        _x = numpy.flip(_profile.y, 0)
        _y = numpy.flip(_profile.x, 0)
        _inverse_profile = xInterpolatedUnivariateSpline(_x, _y, ext='const')

        def _loss(energy, path_length):
            """
            """
            d = _inverse_profile(energy)
            return _profile(d) - _profile(d + path_length)

        return _loss



class xDMEAlphaStoppingPowerTable(xAlphaStoppingPowerTable):

    """Specialized subclass for the DME compound.
    """

    def __init__(self, temperature, pressure):
        """Overloaded constructor.
        """
        self.identifier = 'DME'
        self.density = dme_density(temperature, pressure)
        h = load_alpha_stopping_power_data('H')
        c = load_alpha_stopping_power_data('C')
        o = load_alpha_stopping_power_data('O')
        self._energy = h._energy
        self._stopping_power = h._stopping_power * (6. * H_MASS / DME_MASS) +\
                               c._stopping_power * (2. * C_MASS / DME_MASS) +\
                               o._stopping_power * (1. * O_MASS / DME_MASS)
        self.__calculate_csda_range()
        self.process_data()

    def __calculate_csda_range(self):
        """Integrate the inverse of the stopping power, yielding the range as a
        function of energy in CSDA approximation.
        """
        s = xInterpolatedUnivariateSpline(self._energy, 1. / self._stopping_power)
        self._csda_range = numpy.zeros(self._energy.shape, 'd')
        emin = self._energy.min()
        for i, ebar in enumerate(self._energy):
            self._csda_range[i] = s.integral(emin, ebar)
        self._csda_range += emin * s(0.)



def load_alpha_stopping_power_data(identifier, density=None):
    """Load the stopping power data for a given element or compound.
    """
    return xAlphaStoppingPowerTable(identifier, density)


def load_dme_alpha_stopping_power_data(temperature=20., pressure=800.):
    """Load the stopping power data for a given element or compound.
    """
    return xDMEAlphaStoppingPowerTable(temperature, pressure)

#!/usr/bin/env python
#
# Copyright (C) 2018--2020, the ixpeobssim team.
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


# pylint: disable=invalid-name


class xModelStokesParameters:

    """Small utility class to deal with the Stokes parameters in a source model
    context.

    Basically we provide conversion functions from Stokes parameters to
    polarizarion degree and angle and vice-versa.

    Note that all the algebra, here, is coded in terms of the reduced
    Stokes parameters q = Q / I and u = U / I. The reason is twofold:

    * when we simulate a model we typically decouple the definition of the
      spectrum from that of the polarization pattern;
    * when converting polarization degree and angle to Stokes parameters
      we can only calculate, by definition, q and u---not Q and U.

    For completeness: be aware that all the angles are measured in radians,
    and if you want to operate with degrees it is the user's responsibility
    to do the conversion outside this class.
    """

    @staticmethod
    def polarization_degree(q, u):
        """Convert q and u to the corresponding polarization degree.
        """
        return numpy.sqrt(q**2. + u**2.)

    @staticmethod
    def polarization_angle(q, u):
        """Convert q and u to the corresponding polarization angle (in radians).
        """
        return 0.5 * numpy.arctan2(u, q)

    @staticmethod
    def q(polarization_degree, polarization_angle):
        """Convert a polarization degree and angle (in radians) into the
        corresponding q reduced Stokes parameter.
        """
        return polarization_degree * numpy.cos(2. * polarization_angle)

    @staticmethod
    def u(polarization_degree, polarization_angle):
        """Convert a polarization degree and angle (in radians) into the
        corresponding u reduced Stokes parameter.
        """
        return polarization_degree * numpy.sin(2. * polarization_angle)

    @staticmethod
    def normalize(QU, I):
        """Calcualate the Stokes normalized Q parameter.
        """
        qu = numpy.full(I.shape, 0.)
        mask = I > 0.
        qu[mask] = QU[mask] / I[mask]
        return qu

    @staticmethod
    def pdpa_to_xy(pol_deg, pol_ang, degrees=False):
        """Convert polarization degree and angle into the x and y components
        of the polarization vector.

        This is assuming that the position angle is measured starting from the
        celestial North, see https://bitbucket.org/ixpesw/ixpeobssim/issues/597
        for more discussion about this.

        .. warning::

           This nd the following function are encapsulating our convention for
           measuring position angles, and should be probably better suited in a
           different module?
        """
        if degrees:
            pol_ang = numpy.radians(pol_ang)
        x = -pol_deg * numpy.sin(pol_ang)
        y = pol_deg * numpy.cos(pol_ang)
        return x, y

    @staticmethod
    def qu_to_xy(q, u):
        """Conver the Stokes parameters into the x and y components of the
        polarization vector.
        """
        pd = xModelStokesParameters.polarization_degree(q, u)
        pa = xModelStokesParameters.polarization_angle(q, u)
        return xModelStokesParameters.pdpa_to_xy(pd, pa)



class xDataStokesParameters:

    """Small utility class to deal with the Stokes parameters in a data analysis
    context.

    Warning
    -------
    This is work in progress, and the class methods are known to fail by
    ZeroDivisionError in some circumstances.
    """

    @staticmethod
    def polarization_degree(I, Q, U, dI, dQ, dU):
        """Return the polarization degree and proparate the uncertainties.
        """
        p = xModelStokesParameters.polarization_degree(Q / I, U / I)
        dp = (1. / I) * numpy.sqrt(
            p**2. * dI**2. +
            Q**2. / (I**2. * p**2.) * dQ**2. +
            U**2. / (I**2. * p**2.) * dU**2.)
        return p, dp

    @staticmethod
    def polarization_angle(Q, U, dQ, dU):
        """Return the polarization angle and proparate the uncertainties.
        """
        phi = xModelStokesParameters.polarization_angle(Q, U)
        dphi = numpy.sqrt(
            U**2. / (4 * (Q**2. + U**2)**2.) * dQ**2. +
            Q**2. / (4 * (Q**2. + U**2)**2.) * dU**2.
        )
        return phi, dphi

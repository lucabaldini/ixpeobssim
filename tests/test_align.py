#!/usr/bin/env python
#
# Copyright (C) 2022, the ixpeobssim team.
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

import unittest
import sys

import numpy

from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.core.stokes import xModelStokesParameters
from ixpeobssim.evt.align import align_phi, align_stokes_parameters
from ixpeobssim.evt.kislat2015 import xStokesAnalysis
from ixpeobssim.irf.modf import xAzimuthalResponseGenerator
from ixpeobssim.utils.math_ import fold_angle_rad
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()


class TestAlign(unittest.TestCase):

    """Unit test for the evt.align module in irfgen.
    """

    def test_stokes(self, m=0.5, phase=0., num_events=100000, phi0=0.5):
        """
        """
        binning = numpy.linspace(-numpy.pi, numpy.pi, 100)
        q0 = xStokesAnalysis.stokes_q(phi0, weights=None)
        u0 = xStokesAnalysis.stokes_u(phi0, weights=None)

        # Extract a bunch of randomly distributed phi angles, and calculate
        # the corresponding Stokes parameters.
        m = numpy.full(num_events, m)
        phase = numpy.full(num_events, phase)
        phi = xAzimuthalResponseGenerator().rvs_phi(m, phase)
        q = xStokesAnalysis.stokes_q(phi, weights=None)
        u = xStokesAnalysis.stokes_u(phi, weights=None)
        plt.figure('Original angle')
        xHistogram1d(binning).fill(phi).plot()

        # Plain angle rotation, i.e., old-style phi alignment.
        phi1, q1, u1 = align_phi(phi, phi0)
        plt.figure('Aligned angle')
        xHistogram1d(binning).fill(phi1).plot()

        # Stokes parameter rotation, i.e., new-style Stokes alignment.
        q2, u2 = align_stokes_parameters(q, u, q0, u0)
        # Note that when we do this, phi2 is bound between -pi/2 and pi/2.
        phi2 = xModelStokesParameters.polarization_angle(q2, u2)

        # And the Stokes parameters in the two flavors should be identical.
        self.assertTrue(numpy.allclose(q1, q2))
        self.assertTrue(numpy.allclose(u1, u2))
        self.assertTrue(numpy.allclose(fold_angle_rad(phi1), phi2))



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

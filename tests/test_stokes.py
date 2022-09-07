#!/urs/bin/env python
#
# Copyright (C) 2019, the ixpeobssim team.
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

import numpy

from ixpeobssim.core.stokes import xModelStokesParameters, xDataStokesParameters
from ixpeobssim.utils.logging_ import logger


class TestStokes(unittest.TestCase):

    """Unit test for the utils.stokes classes.
    """

    def base_model_test(self, q, u, polarization_degree, polarization_angle):
        """
        """
        logger.info('Testing Stokes conversion (%.3f, %.3f) <-> (%.3f, %.3f)' %\
                    (q, u, polarization_degree, polarization_angle))
        _pd = xModelStokesParameters.polarization_degree(q, u)
        _pa = xModelStokesParameters.polarization_angle(q, u)
        self.assertAlmostEqual(_pd, polarization_degree)
        self.assertAlmostEqual(_pa, polarization_angle)
        _q = xModelStokesParameters.q(polarization_degree, polarization_angle)
        _u = xModelStokesParameters.u(polarization_degree, polarization_angle)
        self.assertAlmostEqual(_q, q)
        self.assertAlmostEqual(_u, u)

    def test_model(self):
        """
        """
        self.base_model_test(0., 0., 0., 0.)
        self.base_model_test(0.5, 0.5, 1./numpy.sqrt(2.), numpy.pi / 8.)
        self.base_model_test(0.5, 0., 0.5, 0.)

    def test_data(self):
        """
        """
        I, Q, U = 1., 0.1, 0.1
        dI, dQ, dU = 0.01, 0.01, 0.01
        _pd = xDataStokesParameters.polarization_degree(I, Q, U, dI, dQ, dU)
        _pa = xDataStokesParameters.polarization_angle(Q, U, dQ, dU)
        print(_pd, _pa)

    def test_xy(self, pd=1.):
        """Test the orthogonal decomposition of the polarization vector in the sky.

        This was added as a consequence of
        https://bitbucket.org/ixpesw/ixpeobssim/issues/597/
        """
        for pa, targetx, targety in (
            (0., 0., 1.),
            (-0.5 * numpy.pi, 1., 0.),
            (0.5 * numpy.pi, -1., 0.),
            (numpy.pi, 0., -1.)
        ):
            logger.info('Testing decomposition (pd=%.3f, pa=%.3f) -> (x=%.3f, y=%.3f)',\
                pd, pa, targetx, targety)
            x, y = xModelStokesParameters.pdpa_to_xy(1., pa)
            self.assertAlmostEqual(x, targetx)
            self.assertAlmostEqual(y, targety)



if __name__ == '__main__':
    unittest.main()

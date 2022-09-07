#!/usr/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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


"""Unit test for the celestial background components.
"""

import unittest
import sys

import numpy

from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.srcmodel.bkg import xRosatPSPCResponseMatrix, xGalacticBkg
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.utils.matplotlib_ import plt, setup_gca

if sys.flags.interactive:
    plt.ion()



class TestCelestialBkg(unittest.TestCase):

    """Unit test for the celestial background.
    """

    @classmethod
    def setUpClass(cls):
        """
        """
        cls.pspc_resp = xRosatPSPCResponseMatrix()

    def test_pspc_aeff(self):
        """Plot the ROSAT PSPC effective area.
        """
        plt.figure('ROSAT PSPC effective area')
        self.pspc_resp.aeff.plot()
        setup_gca(grids=True)

    def test_rosat_r7_rate(self):
        """Integrate a power law with the spectral index of the galactic
        background, multiplied by the ROSAT PSPC effective area, over the
        standard ROSAT R7 channel (1.05--2.04 keV). This gives the number of
        counts rate per unit area.
        """
        emin = self.pspc_resp.R7_EMIN
        emax = self.pspc_resp.R7_EMAX
        pl_index = xGalacticBkg._SPEC_INDEX
        pl = power_law(1., pl_index)
        aeff = self.pspc_resp.aeff
        energy = numpy.linspace(emin, emax, 100)
        flux = xInterpolatedUnivariateSpline(energy, pl(energy) * aeff(energy))
        counts = flux.integral(emin, emax)
        rad = xGalacticBkg._RADIUS_ARCMIN
        area = 2. * numpy.pi * rad**2.
        norm = 1.e-6 / counts * area
        print(norm)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

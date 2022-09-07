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

import unittest
import os

from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim import IXPEOBSSIM_CONFIG
from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.srcmodel.ephemeris import xOrbitalEphemeris
from ixpeobssim.srcmodel.roi import xCelestialModelComponentBase, xPointSource,\
    xPeriodicPointSource, xUniformDisk, xGaussianDisk, xExtendedSource, \
    xBinarySource, xROIModel
from ixpeobssim.srcmodel.ephemeris import xEphemeris
from ixpeobssim.irf import load_irf_set
import matplotlib.pyplot as plt

"""Unit test for the roi classes.
"""

class TestRoi(unittest.TestCase):

    """Unit test for the ROI classes.
    """

    def test_sources(self):
        """Instantiate all the base classes.
        """
        params = dict(photon_spectrum=power_law(1., 2.),
                      polarization_degree=constant(0.1),
                      polarization_angle=constant(20.),
                      column_density=1.0e-22,
                      redshift=0.0012,
                      identifier=0)
        s1 = xCelestialModelComponentBase('Source 1', **params)
        print(s1)
        s2 = xPointSource('Source 2', ra=10., dec=15., **params)
        print(s2)
        ephem = xEphemeris(10., 0.033)
        s3 = xPeriodicPointSource('Source 3', ra=10., dec=15., ephemeris=ephem, **params)
        print(s3)
        s4 = xUniformDisk('Source 4', ra=10., dec=15., radius=0.01, **params)
        print(s4)
        s5 = xGaussianDisk('Source 5', ra=10., dec=15., sigma=0.01, **params)
        print(s5)
        file_path = os.path.join(IXPEOBSSIM_CONFIG, 'fits', 'crab_0p3_10p0_keV.fits')
        s6 = xExtendedSource('Source 6', file_path, **params)
        print(s6)
        file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'parfiles', 'SAXJ1808.4-3658.par')
        ephem = xOrbitalEphemeris.from_file(file_path)
        s7 = xBinarySource('Source 7', ra=10., dec=15., ephemeris=ephem, **params)
        irf_set = load_irf_set()
        parent_roi = xROIModel(10., 15.)
        print(s7)

if __name__ == '__main__':
    unittest.main()

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
import sys

from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.units_ import atm_to_mbar
from ixpeobssim.irfgen.gpd import dme_photoemission_frac_spline, load_ixpesim_ancillary_data
from ixpeobssim.irfgen.constants import perf_gas_density, dme_density,\
    dme_density_perfect, DME_MASS, DME_REF_DENSITY

if sys.flags.interactive:
    plt.ion()


class TestGpd(unittest.TestCase):

    """Unit test for the GPD module in irfgen.
    """

    def test_dme_density_scaling(self):
        """
        """
        # Test the difference between the calculated and tabulated
        # density values for gaseous DME.
        target_ratio = perf_gas_density(DME_MASS, 0., atm_to_mbar(1.)) / DME_REF_DENSITY
        logger.info('Ratio between calculated and tabulated DME density: %f' %\
                    target_ratio)
        delta = abs(target_ratio - 1.)
        self.assertTrue(delta < 0.05)
        # Test the scaling of the actual DME density with temperature and
        # pressure.
        for temperature in [0, 10, 20, 30]:
            for pressure in [600., 800., 1000.]:
                calc_density = perf_gas_density(DME_MASS, temperature, pressure)
                tab_density = dme_density_perfect(temperature, pressure)
                ratio = calc_density / tab_density
                self.assertAlmostEqual(ratio, target_ratio)

    def test_old_density(self):
        """
        """
        density = perf_gas_density(DME_MASS, 0., atm_to_mbar(1.))
        self.assertAlmostEqual(density, 0.002055, 6)

    @unittest.skip("Disengaged for the public release")
    def test_passive_conversions(self):
        """
        """
        plt.figure('Passive conversions')
        for pressure in [800., 750., 700., 650.]:
            s = dme_photoemission_frac_spline(pressure)
            s.plot(label='%d mbar' % pressure)
        setup_gca(grids=True, legend=True)
        energy, trg_eff, win_prob, gem_prob, _ = load_ixpesim_ancillary_data()
        dme_prob = 1. - win_prob - gem_prob
        s = dme_photoemission_frac_spline(800.)
        for e, p in zip(energy, dme_prob):
            self.assertAlmostEqual(p, s(e))



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

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

import unittest

import numpy

from ixpeobssim.utils.logging_ import logger
import ixpeobssim.config.toy_pollin as srcmod
from ixpeobssim.evt.gti import xUberGTIList
from ixpeobssim.irf import load_irf_set
from ixpeobssim.evt.kislat2015 import xStokesAnalysis
from ixpeobssim.evt.mdp import mdp99
from ixpeobssim.srcmodel.polarization import broadband_pol_deg
from ixpeobssim.utils.time_ import string_to_met_utc



class TestKislat(unittest.TestCase):

    """Unit test for the formalism in Kislat et al., 2015.
    """

    @classmethod
    def setUpClass(cls):
        """Setup the test.
        """
        irf_set = load_irf_set()
        cls.modf = irf_set.modf
        cls.aeff = irf_set.aeff
        kwargs = dict(duration=10000., startdate='2022-04-21', roll=0., gti_list=xUberGTIList())
        kwargs['start_met'] = string_to_met_utc(kwargs.get('startdate'), lazy=True)
        kwargs['stop_met'] = kwargs.get('start_met') + kwargs.get('duration')
        event_list = srcmod.src.rvs_event_list(srcmod, irf_set, **kwargs)
        cls.q, cls.u = event_list.stokes_parameters()
        cls.energy = event_list.energy()
        cls.livetime = kwargs['stop_met'] - kwargs['start_met']

    def _check_delta(self, val, target, threshold=1e-3):
        """Small convenience function to test tolerances.
        """
        delta = abs((val - target) / target)
        logger.info('Value = %.5e, target = %.5e, delta = %.5e', val, target, delta)
        self.assertTrue(delta < threshold)

    def test_naive(self, emin=2., emax=8.):
        """Create a naive, non acceptance-corrected, analysis and make sure it's
        doing the same thing as the traditional machinery.
        """
        logger.info('Testing the un-weighted Stokes analysis...')
        analysis = xStokesAnalysis(self.q, self.u, self.energy, self.modf, self.aeff,
            self.livetime, acceptcorr=False)
        mask = analysis._energy_mask(emin, emax)
        mu = self.modf.weighted_average(analysis._energy[mask])
        counts = numpy.count_nonzero(mask)
        I, Q, U = analysis.sum_stokes_parameters(emin, emax)
        self._check_delta(I, counts)
        mdp = mdp99(mu, counts)
        self._check_delta(analysis.calculate_mdp99(mu, I, analysis.W2(mask)), mdp)
        self._check_delta(analysis.effective_mu(emin, emax), mu)

    def test_polarization(self, emin=2., emax=8.):
        """Run an acceptance-corrected analysis.
        """
        analysis = xStokesAnalysis(self.q, self.u, self.energy, self.modf, self.aeff,
            self.livetime, acceptcorr=True)
        pd, sigma_pd, pa, sigma_pa = analysis.polarization(2., 8., degrees=True)
        self.assertTrue(abs((pa - srcmod.pa) / sigma_pa) < 5.)
        mean_pd = broadband_pol_deg(srcmod.spec, srcmod.pol_deg)
        self.assertTrue(abs((pd - mean_pd) / sigma_pd) < 5.)

    def test_table(self):
        """
        """
        analysis = xStokesAnalysis(self.q, self.u, self.energy, self.modf, self.aeff,
            self.livetime, acceptcorr=False)
        ebinning = numpy.linspace(2., 8., 7)
        table = analysis.polarization_table(ebinning)
        print(table)



if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python
#
# Copyright (C) 2018--2022, the ixpeobssim team.
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

from ixpeobssim.binning.polarization import xBinnedCountSpectrum
from ixpeobssim.core.hist import xScatterPlot
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.irf import load_modf, load_rmf, DEFAULT_IRF_NAME
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()

import ixpeobssim.config.toy_point_source as input_model


class TestStokesSpectra(unittest.TestCase):

    """Unit test for Stokes spectra
    """

    @classmethod
    def setUpClass(cls, mc=True):
        """Setup the test.
        """
        # Initialize the pipeline.
        pipeline.setup(model='toy_point_source', overwrite=True)
        # Generate a test event list.
        evt_list = pipeline.xpobssim(duration=10000.)
        # Bin the data to create PHA1 files (in all three flavors)
        I_list = pipeline.xpbin(*evt_list, algorithm='PHA1', mc=mc)
        U_list = pipeline.xpbin(*evt_list, algorithm='PHA1U', mc=mc)
        Q_list = pipeline.xpbin(*evt_list, algorithm='PHA1Q', mc=mc)
        # Read back the files to get the data.
        cls.I_spec = xBinnedCountSpectrum.from_file_list(I_list)
        cls.U_spec = xBinnedCountSpectrum.from_file_list(U_list)
        cls.Q_spec = xBinnedCountSpectrum.from_file_list(Q_list)
        # Mask the bogus values due to poor statistics
        I = cls.I_spec.RATE
        dI = cls.I_spec.STAT_ERR
        U = cls.U_spec.RATE
        Q = cls.Q_spec.RATE
        mask = (I >= (20 * dI)) * ((U**2. + Q**2.) > 0)
        cls.chan = cls.I_spec.CHANNEL[mask]
        cls.I = cls.I_spec.RATE[mask]
        cls.U = cls.U_spec.RATE[mask]
        cls.Q = cls.Q_spec.RATE[mask]
        modf = load_modf(DEFAULT_IRF_NAME, du_id=1)
        rmf = load_rmf(DEFAULT_IRF_NAME, du_id=1)
        energy = rmf.ebounds(cls.chan)
        cls.mu = modf(energy)

    def test_unitariety(self):
        """
        """
        pd = input_model.pol_deg()
        delta = self.I**2. / (self.U**2. + self.Q**2.) / (pd * self.mu)**2.
        self.assertTrue(abs(delta.mean()) < 0.5 * delta.std())



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

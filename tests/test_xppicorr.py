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
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.irf.ebounds import energy_to_channel
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca

if sys.flags.interactive:
    plt.ion()



class TestPiCorr(unittest.TestCase):

    """Unit test for simulation and analysis pipeline(s).
    """

    @classmethod
    def setUpClass(cls):
        """Create the necessary simluated event file.
        """
        pipeline.reset('toy_point_source', overwrite=True)
        cls.file_list = pipeline.xpobssim(duration=10000.)

    def _test_base(self, label, scale=1.1, **kwargs):
        """Quick test to verify the randomization of the PI while correcting.

        This is generating a simple power-law simulation, runnning xppicorr
        on the output files, and verifying that the scaled PI distribution is
        consistent with what one would obtain by *first* scaling the (floating
        point) energy, and *then* converting back to PI.
        """
        orig_file = xEventFile(self.file_list[0])
        file_list = pipeline.xppicorr(*self.file_list, slope=scale, **kwargs)
        corr_file = xEventFile(file_list[0])
        orig_pi = orig_file.pi_data()
        corr_pi = corr_file.pi_data()
        # Direct energy scaling---note you have to access the column directly,
        # as the energy_data() hook is using the PI.
        energy = orig_file.event_data['ENERGY']
        energy *= scale
        energy_pi = energy_to_channel(energy)
        # Fill the proper histograms.
        binning = numpy.linspace(-0.5, 374.5, 376)
        h1 = xHistogram1d(binning).fill(orig_pi)
        h2 = xHistogram1d(binning).fill(corr_pi)
        h3 = xHistogram1d(binning).fill(energy_pi)
        # Compare the two scaled PI distributions and calculate a chisquare.
        c2 = h2.content
        c3 = h3.content
        mask = (c2 + c3) > 50
        delta = (c2 - c3)[mask] / numpy.sqrt(c2 + c3)[mask]
        chisq = (delta**2.).sum()
        dof = len(delta)
        logger.info('Chisquare: %.3f / %d dof', chisq, dof)

        plt.figure('PI scaling %s' % label)
        h1.plot(label='Original')
        h2.plot(label='Corrected (plain xppicorr)')
        h3.plot(label='Corrected (scaling the energy)')
        setup_gca(logy=True, xlabel='Channel', ylabel='Counts', legend=True)

        plt.figure('PI scaling ratio %s' % label)
        plt.plot(c2 / c3)
        setup_gca(xlabel='Channel', ylabel='Correction ratio', grids=True)

    def test_deterministic(self):
        """
        """
        self._test_base('deterministic')

    def test_non_deterministic(self):
        """
        """
        self._test_base('non deterministic', deterministic=False)





if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

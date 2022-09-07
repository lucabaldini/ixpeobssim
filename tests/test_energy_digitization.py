#!/usr/bin/env python
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


"""Unit test for the energy dispersion facility.
"""

import numpy
import unittest
import sys

from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.irf import load_rmf, DEFAULT_IRF_NAME
from ixpeobssim.irf.ebounds import channel_to_energy, energy_to_channel, ENERGY_STEP
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.core.modeling import xConstant
from ixpeobssim.core.fitting import fit_histogram
from ixpeobssim.core.rand import xUnivariateGenerator
import ixpeobssim.core.pipeline as pipeline

if sys.flags.interactive:
    plt.ion()

class TestEnergyDigitization(unittest.TestCase):

    """Unit test for the energy digitization (added after issue #226).
    """

    def setUp(cls):
        """
        """
        cls.edisp = load_rmf(DEFAULT_IRF_NAME, 1)

    def test_print(self):
        """
        """
        energy = numpy.linspace(1., 1.24, 25)
        pha, pi = self.edisp.pha_analysis(energy)
        for _energy, _pha, _pi in zip(energy, pha, pi):
            print('%.22f' % _energy, _pha, _pi)
        ch = numpy.linspace(0, 10, 11, dtype=numpy.int16)
        energy = self.edisp.channel_to_energy(ch)
        for _ch, _energy in zip(ch, energy):
            print(_ch, _energy)

    def base_test(self, energy, size=1000000):
        """
        """
        ebounds = self.edisp.ebounds
        n = ebounds.num_channels()
        binning = numpy.linspace(0, n, n + 1)
        pha, pi = self.edisp.pha_analysis(energy)
        h1 = xHistogram1d(binning).fill(pi)
        # Now scale the energy.
        scale = numpy.ones(size) - numpy.linspace(0, 0.1, size)
        energy *= scale
        pha, pi = self.edisp.pha_analysis(energy)
        h2 = xHistogram1d(binning).fill(pi)
        return h1, h2

    def test_uniform(self, size=1000000):
        """
        """
        ebounds = self.edisp.ebounds
        energy = numpy.random.uniform(ebounds.min(), ebounds.max(), size=size)
        h1, h2 = self.base_test(energy)
        model = xConstant()
        fit_histogram(model, h1)
        plt.figure('Uniform energy digitization')
        h1.plot()
        model.plot()
        model.stat_box()
        h2.plot()

    def test_power_law(self, size=1000000):
        """
        """
        ebounds = self.edisp.ebounds
        x = numpy.linspace(1., 12., 200)
        y = x ** -2.
        rnd = xUnivariateGenerator(x, y)
        energy = rnd.rvs(size)
        h1, h2 = self.base_test(energy,)
        plt.figure('Power-law energy digitization')
        h1.plot()
        h2.plot()
        plt.yscale('log')
        plt.axis([None, None, 10., None])

    def test_dispersion(self, size=100000):
        """
        """
        ebounds = self.edisp.ebounds
        mc_energy = numpy.random.uniform(ebounds.min(), ebounds.max(), size=size)
        energy, pha, pi = self.edisp.convolve_energy(mc_energy)
        corr_energy, corr_pha, corr_pi = self.edisp.convolve_energy(mc_energy * 0.9)
        n = ebounds.num_channels()
        binning = numpy.linspace(0, n, n + 1)
        h1 = xHistogram1d(binning).fill(pi)
        h2 = xHistogram1d(binning).fill(corr_pi)
        plt.figure('Energy dispersion digitization')
        h2.plot()
        h1.plot()

    def test_conversion(self):
        """Unit test added after https://bitbucket.org/ixpesw/ixpeobssim/issues/531

        We basically want to make sure that the function in irfgen.__init___
        (that are used in xEventList to convert the PI into keV) give exactly
        the same answer as the response matrix class.
        """
        channel = numpy.linspace(0, 374.5, 375).astype(numpy.uint16)
        energy1 = self.edisp.channel_to_energy(channel)
        energy2 = channel_to_energy(channel)
        self.assertTrue(numpy.allclose(energy1, energy2), 'delta = %s' % (energy2 - energy1))

    def test_event_list(self):
        """And this is what we *really* care about: being able to read the
        event energy from the pulse invariant.

        Note that there is a finite resolution (a half of the 40 eV channel step)
        that it's achievable in the conversion.
        """
        pipeline.reset('toy_point_source', overwrite=True)
        file_list = pipeline.xpobssim(duration=100.)
        evt_file = xEventFile(file_list[0])
        energy1 = evt_file.energy_data()
        energy2 = evt_file.event_data['ENERGY']
        max_delta = abs(max(energy1 - energy2))
        self.assertTrue(max_delta <= 0.5 * ENERGY_STEP)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

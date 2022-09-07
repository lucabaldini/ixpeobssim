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


"""Unit test for the irf.spm module.
"""

import unittest
import sys

import numpy

from ixpeobssim.core.fitting import fit_histogram
from ixpeobssim.core.hist import xHistogram1d, xHistogram2d
from ixpeobssim.core.modeling import xModulationCurveRad, xModulationCurveDeg
from ixpeobssim.core.stokes import xModelStokesParameters
from ixpeobssim.irf.modf import xAzimuthalResponseGenerator
from ixpeobssim.irf.spm import xSpuriousModulationMap
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()



class TestSpm(unittest.TestCase):

    """Unit test for the spuriousmodulation maps.
    """

    def _random_positions(self, num_events, radius=2., detx0=0., dety0=0.):
        """Generate random detector positions uniformly within a circle of
        given center and radius.
        """
        r = radius * numpy.sqrt(numpy.random.random(size=num_events))
        phi = 2. * numpy.pi * numpy.random.random(size=num_events)
        detx = detx0 + r * numpy.cos(phi)
        dety = dety0 + r * numpy.sin(phi)
        return detx, dety

    def _test_generation_base(self, true_map, source_q=0., source_u=0.,
                              num_events=1000000, radius=2., detx0=0., dety0=0.):
        """Test generating azimuthal angle distributions.
        """
        detx, dety = self._random_positions(num_events, radius, detx0, dety0)
        phi = true_map.rvs_phi(detx, dety, source_q, source_u)
        binning = numpy.linspace(-numpy.pi, numpy.pi, 200)
        plt.figure('Phi original (%.3f, %.3f)' % (source_q, source_u))
        h = xHistogram1d(binning).fill(phi)
        h.plot()
        model = fit_histogram(xModulationCurveRad(), h)
        model.plot()
        model.stat_box()
        return model

    def _test_generation_constant(self, spur_q=0.1, spur_u=0.):
        """
        """
        spur_mod = xModelStokesParameters.polarization_degree(spur_q, spur_u)
        spur_phase = xModelStokesParameters.polarization_angle(spur_q, spur_u)
        map_ = xSpuriousModulationMap.dummy(spur_q, spur_u, sigma=0.)
        model = self._test_generation_base(map_, 0., 0.)
        model = self._test_generation_base(map_, -spur_q, -spur_u)
        model = self._test_generation_base(map_, spur_q, spur_u)
        model = self._test_generation_base(map_, spur_u, spur_q)

    @staticmethod
    def _plot_phi(phi):
        """
        """
        binning = numpy.linspace(-180., 180., 200)
        h = xHistogram1d(binning, xlabel='$\\phi$ [deg]').fill(numpy.degrees(phi))
        h.plot()
        model = fit_histogram(xModulationCurveDeg(), h)
        model.plot()
        model.stat_box()
        return model

    def _test_correction_base(self, true_map, meas_map, source_q=0., source_u=0.,
                              num_events=1000000, radius=2., detx0=0., dety0=0.):
        """
        """
        source_pol_deg = xModelStokesParameters.polarization_degree(source_q, source_u)
        source_pol_ang = xModelStokesParameters.polarization_angle(source_q, source_u)
        logger.info('Source Stokes parameters: %.5f, %.5f', source_q, source_u)
        logger.info('Source polarization degree: %.5f', source_pol_deg)
        logger.info('Source polarization angle: %.3f degrees', numpy.degrees(source_pol_ang))
        true_map.plot_stokes_params(' true')
        meas_map.plot_stokes_params(' measured')
        detx, dety = self._random_positions(num_events, radius, detx0, dety0)
        phi = true_map.rvs_phi(detx, dety, source_q, source_u)
        phi_corr_true = true_map.correct_phi(phi, detx, dety)
        phi_corr_meas = meas_map.correct_phi(phi, detx, dety)
        plt.figure('Phi original')
        self._plot_phi(phi)
        plt.figure('Phi corrected (true map)')
        model = self._plot_phi(phi_corr_true)
        modulation = model.Modulation
        phase = numpy.radians(model.Phase)
        Q = modulation * numpy.cos(2 * phase)
        U = modulation * numpy.sin(2 * phase)
        logger.info('Delta Q = %.5f, Delta U = %.5f', Q - source_q, U - source_u)
        plt.figure('Phi corrected (measured map)')
        model = self._plot_phi(phi_corr_meas)

    def _test_correction_constant(self, spur_q=0.01, spur_u=0., source_q=0.05, source_u=0.02,
                                 spur_sigma=0.1, num_events=10000000, radius=2.,
                                 detx0=0., dety0=0.):
        """
        """
        true_map = xSpuriousModulationMap.dummy(spur_q, spur_u, sigma=0.)
        meas_map = xSpuriousModulationMap.dummy(spur_q, spur_u, spur_sigma)
        self._test_correction_base(true_map, meas_map, source_q, source_u,
            num_events, radius, detx0, dety0)

    def _test_correction_fabio(self, num_events=25000000, radius=2., detx0=0., dety0=0.):
        """Small test on the source setup proposed by Fabio for the SOC
        pipeline meeting on August 9, 2021.
        """
        source_mod = 0.20
        source_phase = 0.92
        spur_mod = 0.05
        spur_phase = numpy.deg2rad(55.)
        spur_sigma = 0.1
        source_q = source_mod * numpy.cos(2. * source_phase)
        source_u = source_mod * numpy.sin(2. * source_phase)
        spur_q = spur_mod * numpy.cos(2. * spur_phase)
        spur_u = spur_mod * numpy.sin(2. * spur_phase)
        true_map = xSpuriousModulationMap.dummy(spur_q, spur_u, sigma=0.)
        meas_map = xSpuriousModulationMap.dummy(spur_q, spur_u, spur_sigma)
        self._test_correction_base(true_map, meas_map, source_q, source_u,
            num_events, radius, detx0, dety0)

    def _test_correction_split(self, spur_q=0.05, spur_u=0.03, source_q=0.05, source_u=0.02,
                              spur_sigma=0.1, num_events=10000000, radius=2.,
                              detx0=0., dety0=0.):
        """
        """
        Q = numpy.full(xSpuriousModulationMap.SHAPE, spur_q)
        Q[150:, :] = spur_u
        U = numpy.full(xSpuriousModulationMap.SHAPE, spur_u)
        U[150:, :] = spur_q
        true_map = xSpuriousModulationMap.dummy(Q, U, sigma=0.)
        meas_map = xSpuriousModulationMap.dummy(Q, U, spur_sigma)
        self._test_correction_base(true_map, meas_map, source_q, source_u,
            num_events, radius, detx0, dety0)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

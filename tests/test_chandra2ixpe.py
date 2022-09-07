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
import sys

import numpy

from ixpeobssim import IXPEOBSSIM_TEST, IXPEOBSSIM_CONFIG_REG
from ixpeobssim.srcmodel.roi import xChandraObservation, xChandraROIModel
from ixpeobssim.irf import load_irf_set, load_arf, load_vign
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.evt.gti import xUberGTIList
from ixpeobssim.utils.time_ import string_to_met_utc
from ixpeobssim.utils.astro import angular_separation, read_ds9
from ixpeobssim.utils.units_ import degrees_to_arcmin
import ixpeobssim.utils.chandra as chandra
import ixpeobssim.core.pipeline as pipeline

if sys.flags.interactive:
    plt.ion()

class TestChandraToIxpe(unittest.TestCase):

    """Unit test for simulation and analysis pipeline(s).
    """

    def test_chandra_irfs(self):
        """Try to load chandra ACIS-I and ACIS-S arf and vignetting, and
        computing the effective area ratio.
        """
        plt.figure('Chandra effective area')
        acis_i = chandra.load_arf('ACIS-I')
        acis_i.plot(label='ACIS-I')
        acis_s = chandra.load_arf('ACIS-S')
        acis_s.plot(label='ACIS-S')
        plt.legend()
        plt.figure('Chandra vignetting')
        vign_c = chandra.load_vign()
        vign_c.plot()
        plt.figure('Effective area ratio')
        aeff = load_arf()
        vign = load_vign()
        ratio = chandra.arf_ratio(aeff, vign, acis_i, vign_c)
        ratio.plot()

    def test_aeff_ratio(self):
        """Compare the effective area between Chandra and IXPE ratio defined as bivariate
        spline with the one calculated using the effective exposure event by event.
        """
        acis_i = chandra.load_arf('ACIS-I')
        vign_c = chandra.load_vign()
        aeff = load_arf()
        vign = load_vign()
        aeff_ratio = chandra.arf_ratio(aeff, vign, acis_i, vign_c)
        evt_path = os.path.join(IXPEOBSSIM_TEST, 'data', 'cena.fits')
        roi = xChandraROIModel(evt_path, acis='I')
        mask = roi.theta_c < 8.5
        energy = roi.energy_c[mask]
        theta = roi.theta_c[mask]
        ixpe_exp = aeff(energy) * vign(energy, theta) * roi.obs_time
        _ratio1 = ixpe_exp / roi.effexp_c[mask]
        _ratio2 = aeff_ratio(roi.energy_c[mask], roi.theta_c[mask])
        _delta =  _ratio1 / _ratio2 - 1
        _delta_max = numpy.absolute(_delta).max()
        self.assertTrue(_delta_max < 1e-2, 'diff. %.9f' % _delta_max)

    def test_legacy_converter(self, duration=20000.):
        """Compare the new simulation flow using the effective exposure event by event with
        the old one using the time scaling and the effective area ratio.
        """
        numpy.random.seed(0)
        acis_i = chandra.load_arf('ACIS-I')
        vign_c = chandra.load_vign()
        aeff = load_arf()
        vign = load_vign()
        aeff_ratio = chandra.arf_ratio(aeff, vign, acis_i, vign_c)
        evt_path = os.path.join(IXPEOBSSIM_TEST, 'data', 'cena.fits')
        roi = xChandraROIModel(evt_path, acis='I')
        roi._reset_mask()
        mc_energy, mc_ra, mc_dec, mc_effexp, mc_theta =\
            roi.energy_c, roi.ra_c, roi.dec_c, roi.effexp_c, roi.theta_c
        # Old way
        # Scale the chandra events according to duration
        scale = numpy.modf(duration / roi.obs_time)
        _energy, _ra, _dec = xChandraObservation._time_scaling(scale, mc_energy, mc_ra, mc_dec)
        # Throw an array of random numbers and accept the events based on the
        # effective area ratio
        rnd_ratio = numpy.random.random(len(_energy))
        separation = angular_separation(_ra, _dec, roi.ra, roi.dec)
        separation = degrees_to_arcmin(separation)
        mask = rnd_ratio < aeff_ratio(_energy, separation)
        old_num_events = mask.sum()
        # New way
        expect_repeat = aeff(mc_energy) * vign(mc_energy, mc_theta) * duration / mc_effexp
        expect_repeat[numpy.logical_not(expect_repeat > 0.)] = 0.
        new_num_events = expect_repeat.sum()
        _delta =  old_num_events / new_num_events - 1
        _delta_max = numpy.absolute(_delta).max()
        self.assertTrue(_delta_max < 1e-2, 'diff. %.9f' % _delta_max)

    def test_chandra_converter(self):
        """Run a quick test converting data from Chandra.

        Warning
        -------
        We should be checking *something* in the output, here.
        """
        pipeline.reset('test_cena', overwrite=True)
        # Run the converter
        pipeline.xpobssim(duration=10000., saa=False, occult=False)

    def test_nevents(self, N=2):
        """Test if, doubling the duration time, the number of converted events doubles.
        """
        numpy.random.seed(0)
        irf_set = load_irf_set()
        kwargs = dict(startdate='2022-04-21', duration=200000., deadtime=0.,
                      roll=0., dithering=False)
        kwargs['start_met'] = string_to_met_utc(kwargs.get('startdate'), lazy=True)
        kwargs['stop_met'] = kwargs.get('start_met') + kwargs.get('duration')
        kwargs['gti_list'] = xUberGTIList()
        evt_path = os.path.join(IXPEOBSSIM_TEST, 'data', 'cena_mod.fits')
        pol_degree = constant(0.)
        pol_angle = constant(0.)
        obs = xChandraObservation('Cen A', pol_degree, pol_angle)
        roi_model = xChandraROIModel(evt_path, acis='I')
        roi_model.add_source(obs)
        event_list_1 = roi_model.rvs_event_list(irf_set, **kwargs)
        kwargs['duration'] *= N
        event_list_2 = roi_model.rvs_event_list(irf_set, **kwargs)
        _ratio = event_list_2.num_events() / float(event_list_1.num_events())
        _delta = _ratio / N - 1
        self.assertTrue(abs(_delta) < 3e-2, 'diff. %.9f' % abs(_delta))

    def test_overlap(self):
        """Test if the converter raises SystemExit exception in presence of
        overlapping regions in the configuration file.
        """
        irf_set = load_irf_set()
        REG_SOURCE_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG_REG, 'cena_jet+core.reg')
        regions = read_ds9(REG_SOURCE_FILE_PATH)
        kwargs = dict(startdate='2022-04-21', duration=100000., deadtime=0., roll=0.)
        kwargs['start_met'] = string_to_met_utc(kwargs.get('startdate'), lazy=True)
        kwargs['stop_met'] = kwargs.get('start_met') + kwargs.get('duration')
        kwargs['gti_list'] = xUberGTIList()
        evt_path = os.path.join(IXPEOBSSIM_TEST, 'data', 'cena.fits')
        roi_model = xChandraROIModel(evt_path, acis='I')
        pol_degree = constant(0.)
        pol_angle = constant(0.)
        obs1 = xChandraObservation('Cen A', pol_degree, pol_angle)
        roi_model.add_source(obs1)
        obs2 = xChandraObservation('Core', pol_degree, pol_angle, regions[1])
        roi_model.add_source(obs2)
        self.assertRaises(SystemExit, roi_model.rvs_event_list, irf_set,
                          **kwargs)

if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

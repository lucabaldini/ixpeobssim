#!/usr/bin/env python
#
# Copyright (C) 2015, the ixpeobssim team.
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


"""Unit test for the irf.vign module.
"""


import numpy
import unittest

from ixpeobssim.instrument import DU_IDS
from ixpeobssim.core.spline import xInterpolatedBivariateSplineLinear
from ixpeobssim.evt.event import xEventList
from ixpeobssim.irfgen.mma import MMA_VIGN_FILE_PATH
from ixpeobssim.irf import load_irf_set, load_vign, DEFAULT_IRF_NAME



class TestIxpeVign(unittest.TestCase):

    """Unit test for the IXPE vignetting.
    """

    def test_ixpe_vign(self):
        """
        """
        data = numpy.loadtxt(MMA_VIGN_FILE_PATH)
        energy = data.T[0, :]
        vign_data = data[:,1:]
        theta = numpy.arange(0., 9., 0.5)
        _v = xInterpolatedBivariateSplineLinear(energy, theta, vign_data)
        energy, theta = numpy.meshgrid(energy, theta)
        for du_id in DU_IDS:
            vign = load_vign(DEFAULT_IRF_NAME, du_id)
            for _e, _t in zip(energy, theta):
                delta = abs((vign(_e, _t) - _v(_e, _t)) / vign(_e, _t))
                self.assertTrue(delta.max() < 2e-2, 'max. diff. %.9f in DU%d' %\
                                (delta.max(), du_id))

    def test_apply_vignetting(self, size=1000000, rad=0.1):
        """Test the vignetting application on a toy event list.

        This was added for https://bitbucket.org/ixpesw/ixpeobssim/issues/423/
        """
        irf_set = load_irf_set()
        ra0, dec0 = 30., 30.
        evt_list = xEventList()
        for name in evt_list.spec_names:
            evt_list[name] = numpy.zeros(size)
        energy = numpy.full(size, 3.)
        pha, pi = irf_set.edisp.pha_analysis(energy)
        ra = numpy.random.uniform(ra0 - rad, ra0 + rad)
        dec = numpy.random.uniform(dec0 - rad, dec0 + rad)
        evt_list.set_mc_energy_columns(energy, pha, pi)
        evt_list.set_mc_sky_position_columns(ra, dec, 0, 0)
        evt_list.apply_vignetting(irf_set.vign, ra0, dec0)
        frac = evt_list.num_events() / size
        self.assertTrue(frac > 0.1 and frac < 0.95)



if __name__ == '__main__':
    unittest.main()

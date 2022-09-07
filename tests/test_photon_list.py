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

"""Unit tests for photon lists.
"""

from __future__ import print_function, division

import os
import unittest

import numpy

from ixpeobssim import IXPEOBSSIM_DATA
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.evt.gti import xSimpleGTIList
from ixpeobssim.evt.ixpesim import xPhotonList
from ixpeobssim.irf import load_irf_set
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant


# pylint: disable=invalid-name


class TestPhotonList(unittest.TestCase):

    """Unit test for misc.astro module.
    """

    def test(self):
        """
        """
        ra, dec = 30., 45.
        pl_norm = 10.
        pl_index = 2.
        spec = power_law(pl_norm, pl_index)
        pd = 0.1
        pa = 30.
        pol_deg = constant(pd)
        pol_ang = constant(numpy.radians(pa))
        src = xPointSource('Point source', ra, dec, spec, pol_deg, pol_ang)
        roi_model = xROIModel(ra, dec, src)
        irf_set = load_irf_set()
        file_path = os.path.join(IXPEOBSSIM_DATA, 'test_photon_list.fits')
        kwargs = dict(start_met=0., duration=10., gti_list=xSimpleGTIList(0., 10.),
            outfile=file_path)
        photon_list = xPhotonList(numpy.linspace(0., 10., 100), 3)
        photon_list.fill(2.7, 0., 0., 0., 0., 0.1, 0.)
        photon_list.write_fits('Test', roi_model, irf_set, **kwargs)

    def test_ms_pulsar(self):
        """
        """
        pipeline.reset('toy_ms_pulsar', overwrite=True)
        pipeline.xpphotonlist(duration=1000.)


if __name__ == '__main__':
    unittest.main()

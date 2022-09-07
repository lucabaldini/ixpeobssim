#!/usr/bin/env python
#
# Copyright (C) 2020--2022, the ixpeobssim team.
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


"""Unit test for the Crab pulsar model.
"""

import numpy
import unittest
import sys

from ixpeobssim.binning.misc import xBinnedPulseProfile
import ixpeobssim.config.toy_periodic_source as input_model
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt

import ixpeobssim.config.toy_periodic_source as toy_periodic_source

if sys.flags.interactive:
    plt.ion()


"""We explictely set the random seed to have reproducible results.
"""
numpy.random.seed(0)


class TestPulseProfile(unittest.TestCase):

    """Unit test for xCountSpectrum.
    """

    def test_light_curve(self):
        """
        """
        pipeline.reset('toy_periodic_source', overwrite=True)
        file_list = pipeline.xpobssim(duration=1000)
        ephem = toy_periodic_source.ephemeris
        file_list = pipeline.xpphase(*file_list, **ephem.dict())
        file_list = pipeline.xpbin(*file_list, algorithm='PP')
        light_curve = xBinnedPulseProfile.from_file_list(file_list)
        phase = light_curve.PHASE
        counts = light_curve.COUNTS
        error = light_curve.ERROR
        norm = counts.mean() / input_model.pl_norm(phase).mean()
        model = norm * input_model.pl_norm(phase)
        chisq = numpy.sum(((counts - model) / error)**2.)
        ndof = len(model)
        num_sigma = (chisq - ndof) / numpy.sqrt(2. * ndof)
        logger.info('Pulse-profile test: chisq = %.2f / %d (%.2f)', chisq, ndof, num_sigma)
        self.assertTrue(abs(num_sigma) < 5.)
        # Addition for issue #179: make sure the error propagation is correct.
        lc_list = [xBinnedPulseProfile(file_path) for file_path in file_list]
        target_error = numpy.sqrt(sum([lc.ERROR**2. for lc in lc_list]))
        self.assertTrue(numpy.allclose(error, target_error))
        plt.figure('Test pulse profile')
        light_curve.plot()
        plt.plot(phase, model)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

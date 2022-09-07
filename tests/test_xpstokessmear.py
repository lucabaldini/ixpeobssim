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

import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.core.fitting import fit_histogram
from ixpeobssim.core.hist import xHistogram1d, xHistogram2d
from ixpeobssim.core.modeling import xGaussian
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca

if sys.flags.interactive:
    plt.ion()


class TestStokesSmear(unittest.TestCase):

    """Unit test for the smearing of the Stokes parameters.
    """

    def test_toy_point_source(self, sigma=0.025):
        """
        """
        pipeline.reset('toy_point_source', overwrite=True)
        file_list = pipeline.xpobssim(duration=1000.)
        file_path_orig = file_list[0]
        file_list = pipeline.xpstokessmear(file_path_orig, innersigma=sigma, outersigma=sigma)
        file_path_smear = file_list[0]
        evt_file_orig = xEventFile(file_path_orig)
        evt_file_smear = xEventFile(file_path_smear)
        # Retrieve the Stokes parameters.
        qorig, uorig = evt_file_orig.stokes_data()
        qsmear, usmear = evt_file_smear.stokes_data()
        dq = qsmear - qorig
        du = usmear - uorig
        norm_orig = qorig**2. + uorig**2.
        norm_smear = qsmear**2. + usmear**2.
        self.assertTrue(numpy.allclose(norm_orig, 4.))
        self.assertFalse(numpy.allclose(norm_smear, 4.))
        plt.figure('Smeared Stokes parameter normalization')
        xHistogram1d(numpy.linspace(3.5, 4.5, 100)).fill(norm_smear).plot()
        setup_gca('Stokes parameter normalization')
        plt.figure('Delta Stokes parameters')
        hq = xHistogram1d(numpy.linspace(-0.15, 0.15, 100)).fill(dq)
        hu = xHistogram1d(numpy.linspace(-0.15, 0.15, 100)).fill(du)
        hq.plot(label='$\\Delta Q$')
        hu.plot(label='$\\Delta U$')
        gq = fit_histogram(xGaussian(), hq)
        gu = fit_histogram(xGaussian(), hu)
        setup_gca('$\\Delta$ Stokes parameters', legend=True)
        dsigmaq = (gq.parameter_value('Sigma') - sigma) / gq.parameter_error('Sigma')
        dsigmau = (gu.parameter_value('Sigma') - sigma) / gu.parameter_error('Sigma')
        logger.info('Recovered delta sigma Q and U: %.3f, %.3f', dsigmaq, dsigmau)
        self.assertTrue(abs(gq.parameter_value('Sigma') - sigma) < 5.e-3)
        self.assertTrue(abs(gu.parameter_value('Sigma') - sigma) < 5.e-3)
        plt.figure('Smearing delta correlation')
        binning = numpy.linspace(-0.15, 0.15, 50)
        h = xHistogram2d(binning, binning).fill(dq, du)
        h.plot()
        setup_gca(xlabel='$\\Delta$Q', ylabel='$\\Delta$U')



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

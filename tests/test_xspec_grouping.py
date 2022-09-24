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
import os
import glob
import sys

import numpy

from ixpeobssim import IXPEOBSSIM_TEST_DATA
from ixpeobssim.instrument import DU_IDS
import ixpeobssim.config.toy_point_source as srcmod
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.os_ import rm
from ixpeobssim.utils.matplotlib_ import plt, setup_gca

if sys.flags.interactive:
    plt.ion()

from ixpeobssim.utils.environment import PYXSPEC_INSTALLED
if PYXSPEC_INSTALLED:
    from ixpeobssim.evt.xspec_ import compare_fit_data



class TestGrouping(unittest.TestCase):

    """Unit test for rebinning in XSPEC.
    """

    @staticmethod
    def _file_list(*bin_algs):
        """Convenience function to grab the files.
        """
        file_list = []
        for du_id in DU_IDS:
            for bin_alg in bin_algs:
                file_name = 'toy_point_source_du%d_%s.fits' % (du_id, bin_alg)
                file_list.append(os.path.join(IXPEOBSSIM_TEST_DATA, file_name))
        return file_list

    @unittest.skip('Fix issue#382')
    def test_full(self):
        """Test rebinning on a spectro-polarimetric fit with XSPEC.
        """
        if not PYXSPEC_INSTALLED:
            return
        target = {'PhoIndex': srcmod.pl_index, 'norm': srcmod.pl_norm,
                  'poldeg': srcmod.pd, 'polang': srcmod.pa}
        target_vals = list(target.values())
        file_list = self._file_list('pha1', 'pha1q', 'pha1u')
        fit_kwargs = dict(model='powerlaw * polconst', plot=False)
        fit_output = pipeline.xpxspec(*file_list, **fit_kwargs)
        par_names = fit_output.par_names
        par_values = [fit_output.par_values]
        par_errors = [fit_output.par_errors]
        rebin = [1]
        for i in numpy.arange(2, 20):
            kwargs = dict(comm='GROUP 0 275 %d' % i, overwrite=True)
            _file_list = pipeline.xpgrppha(*file_list, **kwargs)
            _fit_output = pipeline.xpxspec(*_file_list, **fit_kwargs)
            par_values.append(_fit_output.par_values)
            par_errors.append(_fit_output.par_errors)
            rebin.append(i)
            self.assertTrue(compare_fit_data(_fit_output, target_vals) == 0)
            if sys.flags.interactive and i == 6:
                pipeline.xpxspec(*_file_list, model='powerlaw * polconst', plot=True)
        for file_path in _file_list:
            rm(file_path)
        for i, par_name in enumerate(par_names):
            plt.figure('rebin %s' % par_name)
            y = [item[i] for item in par_values]
            dy = [item[i] for item in par_errors]
            plt.errorbar(rebin, y, dy, fmt='o')
            plt.axhline(target[par_name], ls='dashed')
            setup_gca(xlabel='Rebin factor', ylabel=par_name)

    @unittest.skip('Fix issue#382')
    def test_normalized(self):
        """Test rebinning on a purely polarimetric fit in XSPEC.
        """
        if not PYXSPEC_INSTALLED:
            return
        target = {'poldeg': srcmod.pd, 'polang': srcmod.pa, 'norm': 1.}
        target_vals = list(target.values())
        file_list = self._file_list('pha1qn', 'pha1un')
        fit_kwargs = dict(model='apolconst', plot=False)
        fit_output = pipeline.xpxspec(*file_list, **fit_kwargs)
        par_names = fit_output.par_names
        par_values = [fit_output.par_values]
        par_errors = [fit_output.par_errors]
        rebin = [1]
        for i in numpy.arange(2, 20):
            kwargs = dict(comm='GROUP 0 275 %d' % i, overwrite=True)
            _file_list = pipeline.xpgrppha(*file_list, **kwargs)
            _fit_output = pipeline.xpxspec(*_file_list, **fit_kwargs)
            par_values.append(_fit_output.par_values)
            par_errors.append(_fit_output.par_errors)
            rebin.append(i)
            self.assertTrue(compare_fit_data(_fit_output, target_vals) == 0)
            if sys.flags.interactive and i == 6:
                pipeline.xpxspec(*_file_list, model='apolconst', plot=True)
        for file_path in _file_list:
            rm(file_path)
        for i, par_name in enumerate(par_names):
            plt.figure('normalized rebin %s' % par_name)
            y = [item[i] for item in par_values]
            dy = [item[i] for item in par_errors]
            plt.errorbar(rebin, y, dy, fmt='o')
            plt.axhline(target[par_name], ls='dashed')
            setup_gca(xlabel='Rebin factor', ylabel=par_name)




if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

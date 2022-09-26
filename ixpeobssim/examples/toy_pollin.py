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

import ixpeobssim.core.pipeline as pipeline
import ixpeobssim.config.toy_pollin as input_model
from ixpeobssim.utils.environment import PYXSPEC_INSTALLED
if PYXSPEC_INSTALLED:
    import ixpeobssim.evt.xspec_ as xspec_


DURATION = 50000.


def simulate():
    """Run the simulation and fold the events in phase.
    """
    pipeline.xpobssim(duration=DURATION)


def bin_():
    """Create the pha1* files.
    """
    for algorithm in ['PHA1', 'PHA1Q', 'PHA1U']:
        pipeline.xpbin(*pipeline.file_list(), algorithm=algorithm)


def spectro_polarimetric_fit():
    """Perform a spectro-polarimetric fit in XSPEC.
    """
    if not PYXSPEC_INSTALLED:
        return
    file_list = pipeline.file_list('pha1*')
    fit_output = pipeline.xpxspec(*file_list, model='powerlaw*pollin',
                                  plot=True, figname=pipeline.figure_name(''),
                                  fixpars='6,0.')
    target = (input_model.pl_index, input_model.pl_norm, input_model.pd_intercept,
              input_model.pd_slope, input_model.pa, 0.)
    xspec_.compare_fit_data(fit_output, target)


def run():
    """Run all.
    """
    simulate()
    bin_()
    spectro_polarimetric_fit()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_pollin')

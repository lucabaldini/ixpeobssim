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

import ixpeobssim.core.pipeline as pipeline
import ixpeobssim.config.toy_offaxis as input_model
from ixpeobssim.binning.exposure import xBinnedLivetimeCube
from ixpeobssim.binning.misc import xBinnedMap
from ixpeobssim.utils.environment import PYXSPEC_INSTALLED
if PYXSPEC_INSTALLED:
    import ixpeobssim.evt.xspec_ as xspec_

DURATION = 20000.


def simulate():
    """Run the simulation.
    """
    pipeline.xpobssim(duration=DURATION, occult=True, saa=True)

def bin_():
    """Create the binned products.
    """
    file_list = pipeline.file_list()
    pipeline.xpbin(*file_list, algorithm='CMAP')
    pipeline.xpbin(*file_list, algorithm='PHA1')
    pipeline.xpbin(*file_list, algorithm='LTCUBE')

def spectral_fit():
    """Fit the PHA1 with and without the correction for the exposure and compare the
    resutls.
    """
    if not PYXSPEC_INSTALLED:
        return
    ltcube_list = pipeline.file_list('ltcube')
    cmap_list = pipeline.file_list('cmap')
    new_arfs = pipeline.xpexposure(*ltcube_list, cmapfiles=cmap_list, irftype='arf')
    pha_list = pipeline.file_list('pha1')
    new_phas = pipeline.xpancrkey(*pha_list, arffiles=new_arfs, overwrite=True)
    fit_output_nocorr = pipeline.xpxspec(*pha_list, model='powerlaw', plot=False)
    fit_output_corr = pipeline.xpxspec(*new_phas, model='powerlaw', plot=False)
    target = (input_model.pl_index, input_model.pl_norm)
    print('Without exposure correction:')
    xspec_.compare_fit_data(fit_output_nocorr, target)
    print('With exposure correction:')
    xspec_.compare_fit_data(fit_output_corr, target)


def display():
    """
    """
    cmap_list = pipeline.file_list('cmap')
    map_ = xBinnedMap.from_file_list(cmap_list)
    pipeline.figure('count map')
    map_.plot()
    ltcube_list = pipeline.file_list('ltcube')
    ltcube = xBinnedLivetimeCube.from_file_list(ltcube_list)
    ltcube.plot()

def run():
    """
    """
    simulate()
    bin_()
    spectral_fit()
    display()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_offaxis')

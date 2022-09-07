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

import os

import astropy
from astropy.io import fits
import numpy

from ixpeobssim import IXPEOBSSIM_BENCHMARKS, IXPEOBSSIM_DATA_BENCHMARKS
from ixpeobssim.benchmarks import plot_pull_histogram
from ixpeobssim.binning.fmt import xBinTableHDUPCUBE
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.core.stokes import xModelStokesParameters
import ixpeobssim.config.toy_point_source as srcmod
from ixpeobssim.utils.matplotlib_ import plt


COL_NAMES = [col_name for col_name, *_ in xBinTableHDUPCUBE.DATA_SPECS]
COL_TYPES = [float] * len(COL_NAMES)
SUMMARY_FILE_NAME = 'stokes_errors_summary.fits'
SUMMARY_FILE_PATH = os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, SUMMARY_FILE_NAME)

#INPUT_PARAMS = (srcmod.pl_index, srcmod.pl_norm, srcmod.pd, srcmod.pa)


def _process(file_list, mc=True):
    """Custom processing routine---we only need the polarization cubes, here.

    Note by default we bin in mc energy, which makes the comparison with
    analytical calculations easier.
    """
    kwargs = dict(emin=2., emax=8., ebins=1, mc=mc)
    pipeline.xpbin(*file_list, algorithm='PCUBE', **kwargs)


def generate(size=1000):
    """Generate an ensamble of realizations.

    This is taylored at using the setup described in
    https://bitbucket.org/ixpesw/ixpeobssim/issues/565
    so that we have solid numbers to compare with.
    """
    kwargs = dict(duration=4502., deadtime=0., vignetting=False, dithering=False)
    pipeline.generate_ensamble(size=size, processing_function=_process, **kwargs)


def _post_process_run(seed=0):
    """Single-run post-processing routine.
    """
    values = []
    file_list = pipeline.glob_ensamble(seed, 'pcube')
    cube = xBinnedPolarizationCube.from_file_list(file_list)
    values = [float(cube.__getattr__(col_name)) for col_name in COL_NAMES]
    return values


def post_process(size=1000):
    """Post process the binned products.
    """
    table = astropy.table.Table(names=COL_NAMES, dtype=COL_TYPES)
    for seed in range(size):
        table.add_row(_post_process_run(seed))
    table.write(SUMMARY_FILE_PATH, format='fits', overwrite=True)


def display(file_path=SUMMARY_FILE_PATH):
    """Display the benchmark results.
    """
    with fits.open(file_path) as hdu_list:
        data = hdu_list[1].data
        #data = {col_name: data[col_name] for col_name in COL_NAMES}

    pd = srcmod.pd
    pa = srcmod.pa
    qn0 = xModelStokesParameters.q(pd, numpy.radians(pa))
    un0 = xModelStokesParameters.u(pd, numpy.radians(pa))

    plt.figure('I pulls')
    plot_pull_histogram(data, 'I')
    plt.figure('QN pulls')
    plot_pull_histogram(data, 'QN', qn0)
    plt.figure('UN pulls')
    plot_pull_histogram(data, 'UN', un0)
    plt.figure('PD pulls')
    plot_pull_histogram(data, 'PD', pd)
    plt.figure('PA pulls')
    plot_pull_histogram(data, 'PA', pa)



def run():
    """Show the results saved into the repo---this isn't re-running anything.
    """
    #file_path = os.path.join(IXPEOBSSIM_BENCHMARKS, 'summary', SUMMARY_FILE_NAME)
    #display(file_path)
    #generate()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_point_source')

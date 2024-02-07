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
import numpy

from ixpeobssim import IXPEOBSSIM_DATA_BENCHMARKS
from ixpeobssim.benchmarks import plot_pull_histogram
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
import ixpeobssim.config.toy_point_source as srcmod
from ixpeobssim.core.hist import xHistogram1d
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.utils.matplotlib_ import plt


COL_NAMES = ['PD0', 'PD0_ERR', 'PD1', 'PD1_ERR', 'PD0MC', 'PD0MC_ERR', 'PD1MC', 'PD1MC_ERR']
COL_TYPES = [float] * len(COL_NAMES)
SUMMARY_FILE_NAME = 'gray_filter_summary.fits'
SUMMARY_FILE_PATH = os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, SUMMARY_FILE_NAME)


def _process(file_list, mc=True):
    """Custom processing routine---we only need the polarization cubes, here.

    We produce the polarization cubes in 6 energy bins, 1 keV wide, binned in
    both true and measured energy. The setup is designed to facilitate the
    comparison with https://github.com/lucabaldini/ixpeobssim/issues/714
    """
    kwargs = dict(emin=2., emax=8., ebins=6, ebinalg='LIN')
    pipeline.xpbin(*file_list, algorithm='PCUBE', mc=True, suffix='pcube_mc', **kwargs)
    pipeline.xpbin(*file_list, algorithm='PCUBE', mc=False, **kwargs)

def generate(size=100):
    """Generate an ensamble of realizations.

    With 100000 s and the gray filter we get about 1.5 M events per run per DU.
    """
    kwargs = dict(duration=100000., grayfilter=True, vignetting=False, dithering=False)
    pipeline.generate_ensamble(size=size, processing_function=_process, **kwargs)

def _post_process_run(seed=0):
    """Single-run post-processing routine.
    """
    values = []
    for suffix in ('pcube', 'pcube_mc'):
        file_list = pipeline.glob_ensamble(seed, suffix)
        pcube = xBinnedPolarizationCube.from_file_list(file_list)
        values += [pcube.PD[0], pcube.PD_ERR[0], pcube.PD[1], pcube.PD_ERR[1]]
    return values

def post_process(size=100):
    """Post process the binned products.
    """
    table = astropy.table.Table(names=COL_NAMES, dtype=COL_TYPES)
    for seed in range(size):
        table.add_row(_post_process_run(seed))
    table.write(SUMMARY_FILE_PATH, format='fits', overwrite=True)

def display(file_path=SUMMARY_FILE_PATH):
    """Display the summary file.
    """
    with astropy.io.fits.open(file_path) as hdu_list:
        data = hdu_list[1].data
    pd = srcmod.pd
    plt.figure('Pulls measured energy 2--3 keV')
    plot_pull_histogram(data, 'PD0', srcmod.pd)
    plt.figure('Pulls measured energy 3--4 keV')
    plot_pull_histogram(data, 'PD1', srcmod.pd)
    plt.figure('Pulls true energy 2--3 keV')
    plot_pull_histogram(data, 'PD0MC', srcmod.pd)
    plt.figure('Pulls true energy 3--4 keV')
    plot_pull_histogram(data, 'PD1MC', srcmod.pd)
    plt.figure('Polarization cube ratio')
    h = xHistogram1d(numpy.linspace(1., 1.25, 50)).fill(data['PD0'] / data['PD0MC'])
    h.plot()


def run():
    """Default target.
    """
    #generate()
    #post_process()
    display()


if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_point_source')

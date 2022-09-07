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

import os

import numpy

from ixpeobssim.benchmarks import xBenchmarkTable, load_benchmark_data
from ixpeobssim.benchmarks import plot_pull_histogram, plot_metrics_comparison
import ixpeobssim.core.pipeline as pipeline
import ixpeobssim.config.toy_point_source as srcmod
from ixpeobssim import IXPEOBSSIM_BENCHMARKS, IXPEOBSSIM_DATA_BENCHMARKS
from ixpeobssim.utils.matplotlib_ import plt


COL_NAMES = ['FIT_INDEX', 'FIT_NORM', 'FIT_PD', 'FIT_PA',
             'FITN_PD', 'FITN_PA',
             'PCUBE_PD', 'PCUBE_PA',
             'MC_PCUBE_PD', 'MC_PCUBE_PA',
             ]
INPUT_PARAMS = (srcmod.pl_index, srcmod.pl_norm, srcmod.pd, srcmod.pa)
SUMMARY_FILE_NAME = 'toy_polconst_summary.fits'
SUMMARY_FILE_PATH = os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, SUMMARY_FILE_NAME)


def generate(size=1000, duration=15000.):
    """Generate an ensamble of realizations.
    """
    pipeline.generate_ensamble(size=size, duration=duration)


def _post_process_run(seed=0):
    """Single-run post-processing routine.
    """
    values = []
    # Fit with a full powerlaw * polconst model.
    kwargs = dict(model='powerlaw * polconst', params=INPUT_PARAMS)
    values += pipeline.fit_ensamble_stokes_spectra(seed, **kwargs)
    # Fit the normalized Stokes parameters with the apollin model.
    kwargs = dict(model='apolconst', params=INPUT_PARAMS[2:])
    values += pipeline.fit_ensamble_stokes_spectra(seed, normalized=True, **kwargs)[:-2]
    # Post-procees the modulation and polarization cubes.
    values += pipeline.post_process_ensamble_pcubes(seed)
    values += pipeline.post_process_ensamble_pcubes(seed, 'mc_pcube')
    return values


def post_process(size=1000):
    """Post process the binned products.
    """
    creator = os.path.basename(__file__)
    table = xBenchmarkTable(COL_NAMES, creator)
    for seed in range(size):
        table.add_row(_post_process_run(seed))
    table.writeto(SUMMARY_FILE_PATH)


def display(file_path=SUMMARY_FILE_PATH):
    """Display the benchmark results.
    """
    _, data = load_benchmark_data(file_path)
    # Display the pulls of the fit parameters.
    targets = list(INPUT_PARAMS) + [srcmod.pd, srcmod.pa] * 4
    for col_name, target in zip(COL_NAMES, targets):
        pipeline.figure('pull_%s' % col_name)
        plot_pull_histogram(data, col_name, target)
    # And some metrics comparison.
    pipeline.figure('Polarization degree')
    col_names = ['FIT_PD', 'PCUBE_PD', 'MC_PCUBE_PD']
    labels = ['XSPEC polconst fit', 'Modulation cube', 'Modulation cube (MC)',
              'Stokes spectra', 'Stokes spectra (MC)']
    plot_metrics_comparison(data, col_names, numpy.linspace(0.06, 0.14, 100), labels,
        xlabel='Polarization degree', legend=True, grids=True)
    pipeline.figure('Polarization angle')
    col_names = ['FIT_PA', 'PCUBE_PD', 'MC_PCUBE_PD']
    plot_metrics_comparison(data, col_names, numpy.linspace(22., 38., 100), labels,
        xlabel='Polarization angle [$^\\circ$]', legend=True, grids=True)


def run():
    """Show the results saved into the repo---this isn't re-running anything.
    """
    file_path = os.path.join(IXPEOBSSIM_BENCHMARKS, 'summary', SUMMARY_FILE_NAME)
    display(file_path)



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_point_source')

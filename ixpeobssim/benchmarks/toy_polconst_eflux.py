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
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.core.modeling import xGaussian
from ixpeobssim.core.fitting import fit_histogram
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, xStatBox
from ixpeobssim.utils.math_ import format_value_error
from astropy.io import fits

COL_NAMES = ['I', 'EFLUX','MC_I', 'MC_EFLUX']
SUMMARY_FILE_NAME = 'toy_pt_src_eflux.fits'
SUMMARY_DATA_FILE_PATH = os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, SUMMARY_FILE_NAME)
SUMMARY_FILE_PATH = os.path.join(IXPEOBSSIM_BENCHMARKS, 'summary', SUMMARY_FILE_NAME)

def plot_residuals_histogram(data, col_name, target=None, span_sigma=7.5, num_bins=100):
    """General-purpose method for plotting the pulls of a parameter estimate.                                                                                        
    """
    # Retrieve the values, the errors and the residuals.
    values = data[col_name]
    errors = data[xBenchmarkTable.error_label(col_name)]
    error = numpy.mean(errors)
    if target is None:
        target = values.mean()
    residuals = values - target
    # Create the histogram.                                                                                                                                                
    try:
        center = residuals.mean()
    except:
        center = 0.
    if (col_name == 'EFLUX' or col_name == 'MC_EFLUX'):
        units = 'erg/cm^2/s'
    elif (col_name == 'I' or col_name == 'MC_I'):
        units ='1/s'
    else: 
        units = ''
    binning = numpy.linspace(center - span_sigma*error, center + span_sigma*error, num_bins)
    label = 'Residuals, %s' % units
    hist = xHistogram1d(binning, xlabel=label).fill(residuals)
    hist.plot()
    # Fit the histogram.                                                                                                                                                  
    try:
        model = fit_histogram(xGaussian(), hist)
    except RuntimeError:
        return
    model.set_plotting_range(binning.min(), binning.max())
    model.plot()
    model.stat_box()
    # Add a custom stat box.                                                                                                                                           
    mean = values.mean()
    mean_err = values.std(ddof=1) / numpy.sqrt(len(values))
    bias = (mean - target) / target
    box = xStatBox()
    box.add_entry(col_name)
    box.add_entry('Average', mean, mean_err)
    box.add_entry('Target', target)
    box.add_entry('Bias: %.2f%%' % (100. * bias))
    box.plot()
    setup_gca(grids=True)
    return hist, model

def generate(size=1000, duration=15000.):
    """Generate an ensamble of realizations.
    """
    pipeline.generate_ensamble(size=size, duration=duration)

def _post_process_run(seed=0):
    """Single-run post-processing routine.
    """
    values = []
    values += pipeline.post_process_ensamble_pcubes_i(seed)
    values += pipeline.post_process_ensamble_pcubes_eflux(seed)    
    values += pipeline.post_process_ensamble_pcubes_i(seed, 'mc_pcube')
    values += pipeline.post_process_ensamble_pcubes_eflux(seed, 'mc_pcube')
    return values


def post_process(size=1000):
    """Post process the binned products.
    """
    creator = os.path.basename(__file__)
    table = xBenchmarkTable(COL_NAMES, creator)    
    for seed in range(size):
        table.add_row(_post_process_run(seed))
    table.writeto(SUMMARY_DATA_FILE_PATH)


def display(file_path=SUMMARY_FILE_PATH):
    """Display the benchmark results.
    """
    _, data = load_benchmark_data(file_path)
    # Display the residuals of the fit parameters.
    targets = [srcmod.i_norm_ltime, srcmod.eflux] * 2
    for col_name, target in zip(COL_NAMES, targets):
        pipeline.figure('residuals_%s' % col_name)
        plot_residuals_histogram(data, col_name, target)

def run():
    """Show the results saved into the repo---this isn't re-running anything.
    """
    #generate(size=1000, duration=15000.)
    #post_process(size=1000)
    file_path = SUMMARY_FILE_PATH
    display(file_path)



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_point_source')

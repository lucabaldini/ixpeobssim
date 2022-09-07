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

import ixpeobssim.core.pipeline as pipeline
import ixpeobssim.config.toy_pollin as srcmod
from ixpeobssim.benchmarks import xBenchmarkTable, load_benchmark_data
from ixpeobssim.benchmarks import plot_pull_histogram, plot_metrics_comparison
from ixpeobssim import IXPEOBSSIM_BENCHMARKS, IXPEOBSSIM_DATA_BENCHMARKS
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color
from ixpeobssim.srcmodel.polarization import broadband_pol_deg
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.core.hist import xHistogram1d



COL_NAMES = ['FIT_INDEX', 'FIT_NORM', 'FIT_PD0', 'FIT_PD1', 'FIT_PA0', 'FIT_PA1',
             'FITN_PD0', 'FITN_PD1', 'FITN_PA0', 'FITN_PA1',
             'FIT_PD', 'FIT_PA',
             'PCUBE_PD', 'PCUBE_PA',
             'MC_PCUBE_PD', 'MC_PCUBE_PA']
INPUT_PARAMS = (srcmod.pl_index, srcmod.pl_norm, srcmod.pd_intercept,
                srcmod.pd_slope, srcmod.pa, 0.)
MEAN_POL_DEG = broadband_pol_deg(srcmod.spec, srcmod.pol_deg)
SUMMARY_FILE_NAME = 'toy_pollin_summary.fits'
SUMMARY_FILE_PATH = os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, SUMMARY_FILE_NAME)



def generate(size=1000, duration=25000.):
    """Generate an ensamble of realizations.
    """
    pipeline.generate_ensamble(size=size, duration=duration)

def _post_process_run(seed=205):
    """Single-run post-processing routine.
    """
    values = []
    # Fit with a full powerlaw * pollin model.
    kwargs = dict(model='powerlaw * pollin', params=INPUT_PARAMS, fixpars='6,0.')
    values += pipeline.fit_ensamble_stokes_spectra(seed, **kwargs)
    # Fit the normalized Stokes parameters with the apollin model.
    kwargs = dict(model='apollin', params=INPUT_PARAMS[2:], fixpars='4,0.')
    values += pipeline.fit_ensamble_stokes_spectra(seed, normalized=True, **kwargs)[:-2]
    # Restricted fit with a polconst model---mind we discard the spectral
    # parameters, as we are only interested in the broadband polarization.
    params = [srcmod.pl_index, srcmod.pl_norm, MEAN_POL_DEG, srcmod.pa]
    kwargs = dict(model='powerlaw * polconst', params=params)
    values += pipeline.fit_ensamble_stokes_spectra(seed, **kwargs)[4:]
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
    targets = list(INPUT_PARAMS) + list(INPUT_PARAMS)[2:] + [MEAN_POL_DEG, srcmod.pa] * 5
    assert len(targets) == len(COL_NAMES)
    for col_name, target in zip(COL_NAMES, targets):
        if col_name in ['FIT_PA1', 'FITN_PA1']:
            # The slope of the polarization angle is fixed, so we don't plot it.
            continue
        pipeline.figure('pull_%s' % col_name)
        plot_pull_histogram(data, col_name, target)
    # And some metrics comparison.
    pipeline.figure('Polarization degree')
    col_names = ['FIT_PD', 'PCUBE_PD', 'MC_PCUBE_PD']
    labels = ['XSPEC polconst fit', 'Polarization cube', 'Polarization cube (MC)',
              'Modulation cube', 'Modulation cube (MC)']
    binning = numpy.linspace(0.04, 0.08, 100)
    # Calculate the average polarization for the pollin XSPEC fit.
    pd = []
    for norm, index, pd0, pd1 in zip(data['FIT_NORM'], data['FIT_INDEX'],
                                     data['FIT_PD0'], data['FIT_PD1']):
        spec = power_law(norm, index)
        pol_deg = lambda E: pd0 + E * pd1
        pd.append(broadband_pol_deg(spec, pol_deg))
    pd = numpy.array(pd)
    hist = xHistogram1d(binning).fill(pd)
    hist.plot(label='XSPEC pollin fit')
    # Calculate the average polarization for the pollin XSPEC normalized fit.
    pdn = []
    for norm, index, pd0, pd1 in zip(data['FIT_NORM'], data['FIT_INDEX'],
                                     data['FITN_PD0'], data['FITN_PD1']):
        spec = power_law(norm, index)
        pol_deg = lambda E: pd0 + E * pd1
        pdn.append(broadband_pol_deg(spec, pol_deg))
    pdn = numpy.array(pdn)
    histn = xHistogram1d(binning).fill(pdn)
    histn.plot(label='XSPEC pollin fit norm.')

    plot_metrics_comparison(data, col_names, binning, labels,
                            xlabel='Polarization degree', legend=True, grids=True)
    plt.axvline(MEAN_POL_DEG, ls='dashed', color='orange')
    # Summary of the different metrics.
    pipeline.figure('Polarization degree precision')
    x = []
    y = []
    dy = []
    for i, (col_name, label) in enumerate(zip(col_names, labels)):
        vals = data[col_name]
        x.append(i)
        y.append(numpy.mean(vals))
        dy.append(numpy.std(vals) / numpy.sqrt(len(vals)))
    x.append(i + 1)
    y.insert(0, hist.mean())
    dy.insert(0, hist.rms() / numpy.sqrt(hist.num_entries()))
    x.append(i + 2)
    y.insert(1, histn.mean())
    dy.insert(1, histn.rms() / numpy.sqrt(hist.num_entries()))
    plt.errorbar(x, y, dy, fmt='o')
    plt.xticks(x, ['XSPEC pollin fit', 'XSPEC pollin fit norm'] + labels, rotation=30., ha='right')
    for _x, _y in zip(x, y):
        delta = 100. * (_y / MEAN_POL_DEG - 1.)
        plt.text(_x, _y, ' $%.1f\\%%$' % delta, color=last_line_color())
    plt.axhline(MEAN_POL_DEG, ls='dashed', color='orange')
    plt.gcf().subplots_adjust(bottom=0.2)
    setup_gca(xmin=-1, xmax=7, grids=True, ylabel='Broadband polarization degree')

def run():
    """Show the results saved into the repo---this isn't re-running anything.
    """
    file_path = os.path.join(IXPEOBSSIM_BENCHMARKS, 'summary', SUMMARY_FILE_NAME)
    display(file_path)




if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_pollin')

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

import os
import glob

import numpy

from ixpeobssim import IXPEOBSSIM_DATA_BENCHMARKS
from ixpeobssim.benchmarks import xBenchmarkTable, load_benchmark_data
from ixpeobssim.benchmarks import plot_pull_histogram, plot_metrics_comparison
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
from ixpeobssim.core.hist import xHistogram1d
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.os_ import mv
from ixpeobssim.utils.matplotlib_ import plt, setup_gca


__model_name__ = 'ixpesim_g2_p700mbar'

PL_INDEX = 2.
PL_NORM = 4.6
POL_DEG = 0.2
POL_ANG = 22.
INPUT_PARAMS = (PL_INDEX, PL_NORM, POL_DEG, POL_ANG)
COL_NAMES = ['FIT_INDEX', 'FIT_NORM', 'FIT_PD', 'FIT_PA',
             'W_FIT_INDEX', 'W_FIT_NORM', 'W_FIT_PD', 'W_FIT_PA',
             'PCUBE_PD', 'PCUBE_PA', 'PCUBE_MDP_99', 'PCUBE_FRAC_W', 'PCUBE_MU_EFF',
             'MC_PCUBE_PD', 'MC_PCUBE_PA', 'MC_PCUBE_MDP_99', 'MC_PCUBE_FRAC_W', 'MC_PCUBE_MU_EFF',
             'W_PCUBE_PD', 'W_PCUBE_PA', 'W_PCUBE_MDP_99', 'W_PCUBE_FRAC_W', 'W_PCUBE_MU_EFF',
             'W_MC_PCUBE_PD', 'W_MC_PCUBE_PA', 'W_MC_PCUBE_MDP_99', 'W_MC_PCUBE_FRAC_W', 'W_MC_PCUBE_MU_EFF'
             ]
SUMMARY_FILE_NAME = 'ixpesim_summary.fits'
SUMMARY_FILE_PATH = os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, SUMMARY_FILE_NAME)
IRF_VERSION = 999
PRESSURE = 700.
DU_ID = 1
UNWEIGHTED_IRF_NAME = 'custom_p%dmbar_v%3d' % (PRESSURE, IRF_VERSION)
WEIGHTED_IRF_NAME = 'customalpha075_p%dmbar_v%3d' % (PRESSURE, IRF_VERSION)
BINNING_DICT = {
    'pha1': dict(algorithm='PHA1', irfname=UNWEIGHTED_IRF_NAME),
    'pha1q': dict(algorithm='PHA1Q', irfname=UNWEIGHTED_IRF_NAME),
    'pha1u': dict(algorithm='PHA1U', irfname=UNWEIGHTED_IRF_NAME),
    'weights_pha1': dict(algorithm='PHA1', irfname=WEIGHTED_IRF_NAME, weights=True),
    'weights_pha1q': dict(algorithm='PHA1Q', irfname=WEIGHTED_IRF_NAME, weights=True),
    'weights_pha1u': dict(algorithm='PHA1U', irfname=WEIGHTED_IRF_NAME, weights=True),
    'pcube': dict(algorithm='PCUBE', ebins=1, irfname=UNWEIGHTED_IRF_NAME),
    '4bin_pcube': dict(algorithm='PCUBE', ebins=4, irfname=UNWEIGHTED_IRF_NAME),
    'mc_pcube': dict(algorithm='PCUBE', mc=True, ebins=1, irfname=UNWEIGHTED_IRF_NAME),
    'mc_4bin_pcube': dict(algorithm='PCUBE', mc=True, ebins=4, irfname=UNWEIGHTED_IRF_NAME),
    'weights_pcube': dict(algorithm='PCUBE', ebins=1, irfname=WEIGHTED_IRF_NAME, weights=True),
    'weights_4bin_pcube': dict(algorithm='PCUBE', ebins=4, irfname=WEIGHTED_IRF_NAME, weights=True),
    'weights_mc_pcube': dict(algorithm='PCUBE', mc=True, ebins=1, irfname=WEIGHTED_IRF_NAME, weights=True),
    'weights_mc_4bin_pcube': dict(algorithm='PCUBE', mc=True, ebins=4, irfname=WEIGHTED_IRF_NAME, weights=True)
}




def bin_photon_lists(folder_path='/data/work/ixpe/data/mdp/'):
    """Loop over all the ixpesim photon lists (processed with xpsimfmt.py) and
    create all the necessary binning products, saving them into
    IXPEOBSSIM_DATA_BENCHMARKS.
    """
    logger.info('Compiling input file list from %s...', folder_path)
    pattern = os.path.join(folder_path, 'ixpesim_spec_gamma2.00_nh0.00e+00_700mbar_r*_recon_simfmt.fits')
    file_list = glob.glob(pattern)
    file_list.sort()
    logger.info('Done, %d file(s) found.', len(file_list))
    for i, file_path in enumerate(file_list):
        for suffix, kwargs in BINNING_DICT.items():
            src = pipeline.xpbin(file_path, **kwargs)[0]
            file_name = '%s_rnd%04d_du%d_%s.fits' % (__model_name__, i, DU_ID, suffix)
            dest = os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, file_name)
            mv(src, dest)


def _xspec_fit(seed, suffix='', **kwargs):
    """
    """
    file_name = '%s_rnd%04d_du%d%s' % (pipeline.model(), seed, DU_ID, suffix)
    base_path = os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, file_name)
    file_list = ['%s_%s.fits' % (base_path, suffix) for suffix in ('pha1', 'pha1q', 'pha1u')]
    fit_output = pipeline.xpxspec(*file_list, **kwargs)
    values = []
    for _, value, error, _, _, _ in fit_output:
        values += [value, error]
    return values


def _read_cube(seed, cls, suffix=''):
    """
    """
    file_name = '%s_rnd%04d_du%d_%s.fits'  % (pipeline.model(), seed, DU_ID, suffix)
    file_path = os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, file_name)
    cube = cls.from_file_list([file_path])
    values = [float(val) for val in cube.polarization()]
    values += [float(cube.MDP_99), 0., float(cube.FRAC_W), 0., float(cube.MU), 0.]
    return values


def _post_process_run(seed=0):
    """Single-run post-processing routine.
    """
    values = []
    # Fit with a full powerlaw * polconst model.
    kwargs = dict(model='powerlaw * polconst', params=INPUT_PARAMS, plot=False)
    # Fit the unweighted Stokes spectra in XSPEC.
    values += _xspec_fit(seed, '', **kwargs)
    # Fit the weighted Stokes spectra in XSPEC.
    values += _xspec_fit(seed, '_weights', **kwargs)
    # Process the polarization cubes.
    values += _read_cube(seed, xBinnedPolarizationCube, 'pcube')
    values += _read_cube(seed, xBinnedPolarizationCube, 'mc_pcube')
    values += _read_cube(seed, xBinnedPolarizationCube, 'weights_pcube')
    values += _read_cube(seed, xBinnedPolarizationCube, 'weights_mc_pcube')
    return values


def post_process(size=100):
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
    for col_name in ('FIT_INDEX', 'W_FIT_INDEX'):
        pipeline.figure('pull_%s' % col_name)
        plot_pull_histogram(data, col_name, PL_INDEX)
    for col_name in ('FIT_PD', 'W_FIT_PD', 'MC_PCUBE_PD', 'W_MC_PCUBE_PD'):
        pipeline.figure('pull_%s' % col_name)
        plot_pull_histogram(data, col_name, POL_DEG)
    pipeline.figure('weight_fraction')
    binning = numpy.linspace(0., 1.05, 100)
    for col_name in ('MC_PCUBE_FRAC_W', 'W_MC_PCUBE_FRAC_W'):
        h = xHistogram1d(binning).fill(data[col_name])
        h.plot(label=col_name.replace('_FRAC_W', ''))
    setup_gca(legend=True, xlabel='$f_w = \\frac{N_{eff}}{N}$', ymax=135., grids=True)
    pipeline.figure('effective_mu')
    binning = numpy.linspace(0., 0.5, 100)
    for col_name in ('MC_PCUBE_MU_EFF', 'W_MC_PCUBE_MU_EFF'):
        h = xHistogram1d(binning).fill(data[col_name])
        h.plot(label=col_name.replace('_MU_EFF', ''))
    setup_gca(legend=True, xlabel='$\\mu_{eff}$', ymax=135., grids=True)
    pipeline.figure('mdp_mc_energy')
    binning = numpy.linspace(0., 0.1, 100)
    for col_name in ('MC_PCUBE_MDP_99', 'W_MC_PCUBE_MDP_99'):
        h = xHistogram1d(binning).fill(data[col_name])
        h.plot(label=col_name.replace('_MDP_99', ''))
    setup_gca(legend=True, xlabel='MDP @ 99% C.L.', ymax=135., grids=True)
    pipeline.figure('mdp_rec_energy')
    binning = numpy.linspace(0., 0.1, 100)
    for col_name in ('PCUBE_MDP_99', 'W_PCUBE_MDP_99'):
        h = xHistogram1d(binning).fill(data[col_name])
        h.plot(label=col_name.replace('_MDP_99', ''))
    setup_gca(legend=True, xlabel='MDP @ 99% C.L.', ymax=135., grids=True)
    r = (data['W_MC_PCUBE_MDP_99'] / data['MC_PCUBE_MDP_99']).mean()
    logger.info('MDP ratio for (MC energy) polarization cubes: %.3f' % r)
    r = (data['W_PCUBE_MDP_99'] / data['PCUBE_MDP_99']).mean()
    logger.info('MDP ratio for (reconstructed energy) polarization cubes: %.3f' % r)


def run():
    """
    """
    pass




if __name__ == '__main__':
    pipeline.bootstrap_pipeline('ixpesim_g2_p700mbar')

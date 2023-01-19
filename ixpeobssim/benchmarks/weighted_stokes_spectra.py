# Copyright (C) 2023, the ixpeobssim team.
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

"""Small benchmark program to verify the correctness of the estimated uncertainities
on the Stokes spectra used for the spectro-polarimetric fits in xspec.

Note that we are only testing the unweighted case, as ixpeobssim is not equipped
to simulate track weights. A full test of the weighted code path requires a
full Geant 4 simulation, and is impractical on a large number of independent
realizations.
"""

import os

from astropy.io import fits
import numpy

from ixpeobssim import IXPEOBSSIM_DATA, IXPEOBSSIM_BENCHMARKS, IXPEOBSSIM_DATA_BENCHMARKS
from ixpeobssim.binning.polarization import xBinnedCountSpectrum
import ixpeobssim.config.toy_point_source as srcmod
from ixpeobssim.core.fitting import fit_gaussian_iterative, fit_histogram
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.core.modeling import xGaussian
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.core.pipeline import stokes_spectra_ensamble_processing, glob_ensamble
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.irf.ebounds import energy_to_channel
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color, save_all_figures


INPUT_PARAMS = (srcmod.pl_index, srcmod.pl_norm, srcmod.pd, srcmod.pa)
FIT_PARAMS = ('INDEX', 'INDEX_ERR', 'NORM', 'NORM_ERR', 'PD', 'PD_ERR', 'PA', 'PA_ERR')
SUMMARY_FILE_NAME = 'weighted_stokes_spectra_summary.fits'
SUMMARY_FILE_PATH = os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, SUMMARY_FILE_NAME)

def _process(file_list):
    """Custom processing function.
    """
    for algorithm in ['PHA1', 'PHA1Q', 'PHA1U']:
        pipeline.xpbin(*file_list, algorithm=algorithm)


def generate(size=100, duration=25000.):
    """Generate an ensamble of realizations.
    """
    pipeline.generate_ensamble(size=size, duration=duration, processing_function=_process)


def _load_pulls(alg='pha1', size=100, aggregate=True):
    """Load the spectral pulls from the PHA1* files.

    Note that this is less straightforward than it might naively seems, as the
    pulls needs to be accumulated on a DU by DU basis---since the three DUs
    have different effective areas, one cannot average the rates across DUs
    to calculate the pulls.

    Arguments
    ---------
    size : int
        The number of random seeds to read in.

    aggregate : bool
        If True, aggregare the pulls for the three DUs at the end.
    """
    rate = {}
    stat_err = {}
    for seed in range(size):
        file_list = glob_ensamble(seed, alg)
        for file_path in file_list:
            spec = xBinnedCountSpectrum(file_path)
            du_id = spec.du_id()
            _rate = spec.RATE
            _stat_err = spec.STAT_ERR
            try:
                rate[du_id] = numpy.vstack((rate[du_id], _rate))
                stat_err[du_id] = numpy.vstack((stat_err[du_id], _stat_err))
            except KeyError:
                rate[du_id] = _rate
                stat_err[du_id] = _stat_err
    average_rate = {du_id: rate[du_id].mean(axis=0) for du_id in DU_IDS}
    pulls = {du_id: (rate[du_id] - average_rate[du_id]) / stat_err[du_id] for du_id in DU_IDS}
    if aggregate:
        pulls = numpy.vstack(list(pulls.values()))
    return pulls


def _fit(size=100):
    """Loop over the Stokes spectra and perform a spectro-polarimetric fit
    with the proper model.

    This returns a 2-dimensional array, each row containing a full set of
    fit parameters with the associated errors.
    """
    fit_results = []
    kwargs = dict(model='powerlaw*polconst', params=INPUT_PARAMS)
    for seed in range(size):
        fit_results.append(pipeline.fit_ensamble_stokes_spectra(seed, **kwargs))
    return numpy.array(fit_results)


def process(size=100):
    """Process the Stokes spectra and write the output FITS file with all the
    summary data.
    """
    logger.info('Processing Stokes spectra...')
    ipulls = _load_pulls('pha1', size)
    qpulls = _load_pulls('pha1q', size)
    upulls = _load_pulls('pha1u', size)
    fit_results = _fit(size)
    logger.info('Preparing FITS summary file...')
    pull_fmt = '375E'
    pull_cols = (
        fits.Column(name='IPULLS', format=pull_fmt, array=ipulls),
        fits.Column(name='QPULLS', format=pull_fmt, array=qpulls),
        fits.Column(name='UPULLS', format=pull_fmt, array=upulls)
    )
    pull_hdu = fits.BinTableHDU.from_columns(pull_cols)
    pull_hdu.header.set('EXTNAME', 'PULLS')
    fit_cols = [fits.Column(name=col_name, format='E', array=fit_results[:,i]) for i, col_name in enumerate(FIT_PARAMS)]
    fit_hdu = fits.BinTableHDU.from_columns(fit_cols)
    fit_hdu.header.set('EXTNAME', 'FIT')
    hdu_list = fits.HDUList([fits.PrimaryHDU(), pull_hdu, fit_hdu])
    logger.info('Writing summary file to %s...', SUMMARY_FILE_PATH)
    hdu_list.writeto(SUMMARY_FILE_PATH, overwrite=True)


def _plot_pulls(pull_name, pulls, emin, emax):
    """Small convenience function to plot the pulls.
    """
    chmin = int(energy_to_channel(emin))
    chmax = int(energy_to_channel(emax)) + 1
    pulls = pulls[:,chmin:chmax]
    pull_binning = numpy.linspace(-5., 5., 100)
    plt.figure('%s pull distribution' % pull_name)
    h = xHistogram1d(pull_binning, xlabel='%s pulls' % pull_name).fill(pulls.flatten())
    h.plot()
    model = fit_histogram(xGaussian(), h)
    model.plot()
    model.stat_box()
    plt.figure('%s pull profile' % pull_name)
    ch = numpy.arange(chmin, chmax)
    plt.plot(ch, pulls.mean(axis=0), label='Pull mean')
    plt.axhline(0., color=last_line_color(), ls='dashed')
    plt.plot(ch, pulls.std(axis=0), label='Pull sigma')
    plt.axhline(1., color=last_line_color(), ls='dashed')
    setup_gca(xlabel='Channel', legend=True, xmin=chmin, xmax=chmax)


def _plot_fit_pulls(fit_data):
    """Plot the pulls in the fit parameters.
    """
    pull_binning = numpy.linspace(-10., 10., 100)
    keys = [key for key in FIT_PARAMS if '_ERR' not in key]
    for key, target in zip(keys, INPUT_PARAMS):
        plt.figure('%s fit pulls' % key)
        pulls = (fit_data[key] - target) / fit_data['%s_ERR' % key]
        h = xHistogram1d(pull_binning, xlabel='%s pulls' % key).fill(pulls)
        h.plot()
        model = fit_histogram(xGaussian(), h)
        model.plot()
        model.stat_box()


def display(file_path, emin=2., emax=8.):
    """Plot all the relevant stuff.
    """
    logger.info('Reading summary file from %s...', file_path)
    hdu_list = fits.open(file_path)
    _plot_pulls('I', hdu_list['PULLS'].data['IPULLS'], emin, emax)
    _plot_pulls('Q', hdu_list['PULLS'].data['QPULLS'], emin, emax)
    _plot_pulls('U', hdu_list['PULLS'].data['UPULLS'], emin, emax)
    _plot_fit_pulls(hdu_list['FIT'].data)


def run():
    """Default entry point.
    """
    file_path = os.path.join(IXPEOBSSIM_BENCHMARKS, 'summary', SUMMARY_FILE_NAME)
    display(file_path)
    save_all_figures(IXPEOBSSIM_DATA)
    plt.show()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_point_source')

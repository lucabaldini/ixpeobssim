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


import numpy

from ixpeobssim.binning.polarization import xBinnedCountSpectrum
from ixpeobssim.core.fitting import fit_gaussian_iterative
from ixpeobssim.core.hist import xHistogram1d
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.core.pipeline import stokes_spectra_ensamble_processing, glob_ensamble
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.irf.ebounds import energy_to_channel
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca



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
    average_rate = {du_id: rate[du_id].mean(axis=0) for du_id in rate.keys()}
    pulls = {du_id: (rate[du_id] - average_rate[du_id]) / stat_err[du_id] for du_id in rate.keys()}
    if aggregate:
        pulls = numpy.vstack(list(pulls.values()))
    return pulls


def post_process(size=100, emin=2., emax=8.):
    """Read back the binned count spectra and create the relevant histograms
    of the pulls.
    """
    chmin = int(energy_to_channel(emin))
    chmax = int(energy_to_channel(emax)) + 1
    for alg in ('pha1', 'pha1q', 'pha1u'):
        pulls = _load_pulls(alg, size)
        print(pulls.shape)
        plt.figure('%s Pulls' % alg.upper())
        h = xHistogram1d(numpy.linspace(-5., 5., 100))
        vals = pulls[:,chmin:chmax].flatten()
        h.fill(vals)
        h.plot()
        setup_gca(xlabel='%s Pulls' % alg.upper())
        try:
            model = fit_gaussian_iterative(h, num_sigma_left=3., num_sigma_right=3.)
            model.plot()
            model.stat_box()
        except RuntimeError:
            pass
    plt.show()


def run():
    """Default entry point.
    """
    pass



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_point_source')

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

import numpy

import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.utils.misc import pairwise, pairwise_enum
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, residual_plot
import ixpeobssim.config.toy_livetime_phase as input_model


def simulate():
    """Run the simulation and fold the events in phase.
    """
    kwargs = dict(occult=False, saa=False)
    file_list = pipeline.xpobssim(duration=50000, **kwargs)
    pipeline.xpphase(*file_list, suffix='folded', **input_model.ephemeris.dict())

def _suffix(i):
    """Convenience function to label phase selections.
    """
    return pipeline.suffix('phase', i)

def select():
    """Run xpselect in a number of bins in phase.
    """
    file_list = pipeline.file_list('folded')
    for i, (tmin, tmax) in pairwise_enum(input_model.PHASE_INTERVALS):
        kwargs = dict(phasemin=tmin, phasemax=tmax, suffix=_suffix(i), ltimeupdate=True, ltimealg='LTSCALE')
        pipeline.xpselect(*file_list, **kwargs)

def bin_(rebin=20):
    """Create all the necessary binned products for a spectro-polarimetric fit.
    """
    for algorithm in ('PHA1', 'PHA1Q', 'PHA1U'):
        for i in range(input_model.NUM_INTERVALS):
            file_list = pipeline.file_list('folded', _suffix(i))
            file_list = pipeline.xpbin(*file_list, algorithm=algorithm)
            pipeline.xpgrppha(*file_list, comm='GROUP 0 375 %d' % rebin)

def fit_phase_bin(i=0):
    """Run an XSPEC spectro-polarimetric fit for a single phase bin.
    """
    file_list = pipeline.file_list('folded', _suffix(i), 'pha1*', 'grppha')
    fit_results = pipeline.xpxspec(*file_list, model='polconst*powerlaw', plot=False)
    return fit_results

def display():
    """Fit all the phase bins and display the fit parameters as a function of the
    phase bin.
    """
    fit_results = [fit_phase_bin(i) for i in range(input_model.NUM_INTERVALS)]
    phase = numpy.array([0.5 * (phasemin + phasemax) for \
        (phasemin, phasemax) in pairwise(input_model.PHASE_INTERVALS)])
    phase_grid = numpy.linspace(0., 1., 1000)
    index = lambda p: numpy.full(p.shape, input_model.PL_INDEX)
    pol_deg = lambda p:numpy.full(p.shape, input_model.POL_DEG)
    pol_ang = lambda p: numpy.full(p.shape, input_model.POL_ANG)
    input_ = index, input_model.pl_norm, pol_deg, pol_ang
    for par_name, model in zip(('PhoIndex', 'norm', 'A', 'psi'), input_):
        val = numpy.array([res.par_value(par_name) for res in fit_results])
        err = numpy.array([res.par_error(par_name) for res in fit_results])
        ax1, ax2 = residual_plot(par_name)
        plt.errorbar(phase, val, err, fmt='o', label='XSPEC fit')
        plt.plot(phase_grid, model(phase_grid), label='Input model')
        setup_gca(ylabel=par_name, grids=True, legend=True)
        plt.sca(ax2)
        plt.errorbar(phase, val - model(phase), err, fmt='o')
        setup_gca(xmin=0., xmax=1., xlabel='Phase', ylabel=par_name, grids=True)
    plt.show()

def run():
    """Run all.
    """
    simulate()
    select()
    bin_()
    display()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_livetime_phase')

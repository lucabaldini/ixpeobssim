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
import ixpeobssim.config.toy_livetime_time as input_model


def simulate():
    """Run the simulation.
    """
    kwargs = dict(startdate=input_model.START_DATE, duration=input_model.DURATION,
        occult=False, saa=False)
    pipeline.xpobssim(**kwargs)

def _suffix(i):
    """Convenience function to label time selections.
    """
    return pipeline.suffix('interval', i)

def select():
    """Run xpselect in a number of bins in time.
    """
    file_list = pipeline.file_list()
    for i, (tmin, tmax) in pairwise_enum(input_model.TIME_INTERVALS):
        kwargs = dict(tmin=tmin, tmax=tmax, suffix=_suffix(i), ltimeupdate=True, ltimealg='LTSUM')
        pipeline.xpselect(*file_list, **kwargs)

def bin_(rebin=20):
    """Create all the necessary binned products for a spectro-polarimetric fit.
    """
    for algorithm in ('PHA1', 'PHA1Q', 'PHA1U'):
        for i in range(input_model.NUM_INTERVALS):
            file_list = pipeline.file_list(_suffix(i))
            file_list = pipeline.xpbin(*file_list, algorithm=algorithm)
            pipeline.xpgrppha(*file_list, comm='GROUP 0 375 %d' % rebin)

def fit_time_interval(i=0):
    """Run an XSPEC spectro-polarimetric fit for a single time interval.
    """
    file_list = pipeline.file_list(_suffix(i), 'pha1*', 'grppha')
    fit_results = pipeline.xpxspec(*file_list, model='polconst*powerlaw', plot=False)
    return fit_results

def display():
    """Fit all the phase bins and display the fit parameters as a function of the
    phase bin.
    """
    fit_results = [fit_time_interval(i) for i in range(input_model.NUM_INTERVALS)]
    time_ = numpy.array([0.5 * (tmin + tmax) for (tmin, tmax) in pairwise(input_model.TIME_INTERVALS)])
    time_grid = numpy.linspace(min(time_), max(time_), 1000)
    index = lambda t: numpy.full(t.shape, input_model.PL_INDEX)
    pol_deg = lambda t:numpy.full(t.shape, input_model.POL_DEG)
    pol_ang = lambda t: numpy.full(t.shape, input_model.POL_ANG)
    input_ = index, input_model.pl_norm, pol_deg, pol_ang
    for par_name, model in zip(('PhoIndex', 'norm', 'A', 'psi'), input_):
        val = numpy.array([res.par_value(par_name) for res in fit_results])
        err = numpy.array([res.par_error(par_name) for res in fit_results])
        ax1, ax2 = residual_plot(par_name)
        plt.errorbar(time_, val, err, fmt='o', label='XSPEC fit')
        plt.plot(time_grid, model(time_grid), label='Input model')
        setup_gca(ylabel=par_name, grids=True, legend=True)
        plt.sca(ax2)
        plt.errorbar(time_, val - model(time_), err, fmt='o')
        setup_gca(xlabel='time_', ylabel=par_name, grids=True)
    plt.show()

def run():
    """Run all.
    """
    simulate()
    select()
    bin_()
    display()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_livetime_time')

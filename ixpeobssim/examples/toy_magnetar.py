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
from ixpeobssim.binning.base import xEventBinningBase
from ixpeobssim.binning.misc import xBinnedPulseProfile
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
from ixpeobssim.utils.misc import pairwise_enum
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color
from ixpeobssim.utils.fmtaxis import fmtaxis, label
import ixpeobssim.config.toy_magnetar as input_model



DURATION = 200000.
PHASE_BINNING = numpy.linspace(0., 1., 25)
ENERGY_BINNING = numpy.array([2., 8.])
USE_MC = True


def simulate():
    """Run the simulation and fold the events in phase.
    """
    file_list = pipeline.xpobssim(duration=DURATION)
    pipeline.xpphase(*file_list, suffix='folded', **input_model.ephemeris.dict())

def select():
    """Run xpselect in a number of bins in phase.
    """
    file_list = pipeline.file_list('folded')
    for i, (min_, max_) in pairwise_enum(PHASE_BINNING):
        pipeline.xpselect(*file_list, phasemin=min_, phasemax=max_,
                          suffix=pipeline.suffix('phase', i))

def bin_():
    """Create a pulse profile, as well as modulation a modulation for each
    subselection in phase.
    """
    # Pulse profile.
    file_list = pipeline.file_list('folded')
    pipeline.xpbin(*file_list, algorithm='PP')
    for i, (min_, max_) in pairwise_enum(PHASE_BINNING):
        file_list = pipeline.file_list('folded', ('phase', i))
        pipeline.xpbin(*file_list, algorithm='PCUBE', ebinalg='LIST',
                       ebinning=ENERGY_BINNING, mc=USE_MC)

def display_pulse_profile():
    """Display the pulse profile.
    """
    pipeline.figure('pulse profile')
    file_list = pipeline.file_list('folded_pp')
    pp = xBinnedPulseProfile.from_file_list(file_list)
    pp.plot(label='IXPE %d ks' % pp.ontime())
    phase = numpy.linspace(0., 1., 100)
    model_pp = numpy.full(phase.shape, numpy.mean(pp.COUNTS))
    plt.plot(phase, model_pp, label='Input model')
    setup_gca(legend=True)


def display_polarization():
    """Display the polarization degree as a function of the pulse phase.
    """
    phase = numpy.linspace(0., 1., 100)
    phase_bins = xEventBinningBase.bin_centers(PHASE_BINNING)
    shape = (len(ENERGY_BINNING) - 1, len(PHASE_BINNING) - 1)
    pol_deg = numpy.zeros(shape)
    pol_deg_err = numpy.zeros(shape)
    pol_ang = numpy.zeros(shape)
    pol_ang_err = numpy.zeros(shape)
    emean = numpy.zeros(shape)
    for i, (min_, max_) in pairwise_enum(PHASE_BINNING):
        file_list = pipeline.file_list('folded', ('phase', i), 'pcube')
        pcube = xBinnedPolarizationCube.from_file_list(file_list)
        pol_deg[:,i] = pcube.PD
        pol_deg_err[:,i] = pcube.PD_ERR
        pol_ang[:,i] = pcube.PA
        pol_ang_err[:,i] = pcube.PA_ERR
        emean[:,i] = pcube.E_MEAN

    def data_label(emin, emax):
        return 'IXPE %d ks (%.2f - %.2f keV)' % (pcube.ontime(), emin, emax)

    def model_label():
        return 'Input model'

    ax1, ax2 = pipeline.residual_figure('polarization degree')
    for i, (min_, max_) in pairwise_enum(ENERGY_BINNING):
        pd = pol_deg[i,:]
        pderr = pol_deg_err[i,:]
        modpd = input_model.pol_deg(None, phase)
        plt.sca(ax1)
        plt.errorbar(phase_bins, pd, pderr, fmt='o', label=data_label(min_, max_))
        llc = last_line_color()
        plt.sca(ax2)
        plt.errorbar(phase_bins, pd - input_model.pol_deg(None, phase_bins), pderr,
                     fmt='o', color=llc)
        plt.sca(ax1)
        plt.plot(phase, modpd, color=llc, label=model_label())
    plt.sca(ax1)
    setup_gca(ymin=0, ymax=1., legend=True, **fmtaxis.pp_pol_deg)
    plt.sca(ax2)
    plt.hlines(0., 0., 1.)
    ymin, ymax = plt.ylim()
    dy = max(abs(ymin), abs(ymax))
    setup_gca(ymin=-dy, ymax=dy, xlabel=label.phase, ylabel='Residuals')

    pipeline.figure('polarization angle')
    for i, (min_, max_) in pairwise_enum(ENERGY_BINNING):
        plt.errorbar(phase_bins, pol_ang[i,:], pol_ang_err[i,:], fmt='o',
                     label=data_label(min_, max_))
        energy = numpy.mean(emean[i,:])
        plt.plot(phase, numpy.degrees(input_model.pol_ang(energy, phase)),
                 color=last_line_color(), label=model_label())
    setup_gca(legend=True, xlabel=label.phase)


def display():
    """Display all.
    """
    display_pulse_profile()
    display_polarization()


def run():
    """Run all.
    """
    simulate()
    select()
    bin_()
    display()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_magnetar')

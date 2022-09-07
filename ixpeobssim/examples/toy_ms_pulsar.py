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
import ixpeobssim.config.toy_ms_pulsar as input_model


# Common settings for the simulation.
DURATION = 25000.
SIM_KWARGS = dict(occult=False, saa=False, duration=DURATION, seed=101)


def create_event_list():
    """Run the simulation and fold the events in phase.
    """
    file_list = pipeline.xpobssim(**SIM_KWARGS)
    pipeline.xpphase(*file_list, suffix='folded', **input_model.ephemeris.dict())

def create_photon_list():
    """Create a photon list to be fed into ixpesim. This can be fed into ixpesim
    with something along the lines of

    ixpesim
        --src-photon-list ~/ixpeobssimdata/toy_ms_pulsar_du1_photon_list.fits
        -n 100000000
        --dme-pressure 645.
        --output-file ~/ixpeobssimdata/toy_ms_pulsar_du1_ixpesim.fits
        --dead-time-offset 760
        --dead-time-slope 300

    See https://bitbucket.org/ixpesw/gpdsw/issues/334 for all the indication
    about the deadtime.
    """
    file_list = pipeline.xpphotonlist(**SIM_KWARGS)

def run():
    """Run all.
    """
    create_event_list()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_ms_pulsar')

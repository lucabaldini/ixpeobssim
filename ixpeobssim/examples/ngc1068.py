# Copyright (C) 2018-2022, the ixpeobssim team.
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
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
from ixpeobssim.utils.matplotlib_ import setup_gca

import ixpeobssim.config.ngc1068 as input_model


DURATION = 2500000.
ENERGY_BINNING = numpy.array([2., 2.75, 4., 6., 8.])


def simulate():
    """Run the simulation.
    """
    pipeline.xpobssim(duration=DURATION)


def bin_():
    """Create the binned pcube object.
    """
    file_list = pipeline.file_list()
    pipeline.xpbin(*file_list, algorithm='PCUBE', ebinalg='LIST',
                   ebinning=ENERGY_BINNING)


def display():
    """Display the stuff.
    """
    file_list = pipeline.file_list('pcube')
    pcube = xBinnedPolarizationCube.from_file_list(file_list)
    ontime = 1.e-3 * pcube.ontime()
    pipeline.figure('polarization degree')
    pcube.plot_polarization_degree(label='IXPE %.1f Ms' % ontime)
    input_model.pol_deg_spline.plot(label='Input model')
    setup_gca(xmin=1., xmax=10., ymin=0., legend=True)


def run():
    """Run the pipeline.
    """
    simulate()
    bin_()
    display()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('ngc1068')

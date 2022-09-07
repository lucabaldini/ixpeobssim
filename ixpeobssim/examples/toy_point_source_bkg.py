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
import ixpeobssim.config.toy_point_source_bkg as input_model
from ixpeobssim.binning.polarization import xBinnedPolarizationCube


DURATION = 5000000.
ENERGY_BINNING = numpy.array([2., 3., 4., 5.5, 8.])


def simulate():
    """Run the simulation.
    """
    pipeline.xpobssim(duration=DURATION, saa=False, occult=False)


def select(src_rad=0.75, bkg_inner_rad=1.5, bkg_outer_rad=3.):
    """Select the photon lists.
    """
    file_list = pipeline.file_list()
    pipeline.xpselect(*file_list, rad=src_rad, suffix='src')
    pipeline.xpselect(*file_list, innerrad=bkg_inner_rad, rad=bkg_outer_rad, suffix='bkg')


def bin_():
    """Create the necessary binned files.
    """
    pipeline.xpbin(*pipeline.file_list(), algorithm='CMAP')
    pipeline.xpbin(*pipeline.file_list('src'), algorithm='CMAP')
    pipeline.xpbin(*pipeline.file_list('bkg'), algorithm='CMAP')
    kwargs = dict(algorithm='PCUBE', ebinalg='LIST', ebinning=ENERGY_BINNING)
    pipeline.xpbin(*pipeline.file_list('src'), **kwargs)
    pipeline.xpbin(*pipeline.file_list('bkg'), **kwargs)

def subtract():
    """
    """
    pcube_src = xBinnedPolarizationCube.from_file_list(pipeline.file_list('src', 'pcube'))
    pd_original = pcube_src.PD
    pd_err_original = pcube_src.PD_ERR
    pcube_bkg = xBinnedPolarizationCube.from_file_list(pipeline.file_list('bkg', 'pcube'))
    beta = pcube_src.backscal() / pcube_bkg.backscal()
    pcube_bkg *= beta
    pcube_src -= pcube_bkg
    print(pcube_src.PD)
    print(pcube_src.PD_ERR)
    print((pcube_src.PD - input_model.PD) / pcube_src.PD_ERR)
    print(pd_original)
    print(pd_err_original)

def run():
    """Run all.
    """
    simulate()
    select()
    bin_()


if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_point_source_bkg')

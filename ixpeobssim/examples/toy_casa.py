# Copyright (C) 2016--2022, the ixpeobssim team.
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

"""Analysis pipeline for the toy_casa example.
"""

from __future__ import print_function, division

import os

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.binning.misc import xBinnedMap
from ixpeobssim.binning.polarization import xBinnedPolarizationMapCube
from ixpeobssim.core import pipeline


DURATION = 1500000.
ENERGY_BINNING = [2., 4., 8.]


def simulate():
    """Run the simulation and align the stokes parameters tangentially.
    """
    file_list = pipeline.xpobssim(duration=DURATION)
    pipeline.xpstokesalign(*file_list, mode='TAN')

def bin_():
    """Create the counts map and the polarization map cubes, in both the original
    and the rotated (via xpstokealign) flavors.
    """
    kwargs = dict(npix=50, ebinalg='LIST', ebinning=ENERGY_BINNING)
    file_list = pipeline.file_list()
    pipeline.xpbin(*file_list, algorithm='CMAP')
    pipeline.xpbin(*file_list, algorithm='PMAPCUBE', **kwargs)
    file_list = pipeline.file_list('stokesalign')
    pipeline.xpbin(*file_list, algorithm='PMAPCUBE', **kwargs)

def display():
    """Display the binned data products---this will plot the counts map and the
    polarization cubes in the two flavors.
    """
    # Plot the count map.
    pipeline.figure('count map')
    file_list = pipeline.file_list('cmap')
    map_ = xBinnedMap.from_file_list(file_list)
    map_.plot()
    # Plot the polarization map cube before the alignment.
    file_list = pipeline.file_list('pmapcube')
    pmap_cube = xBinnedPolarizationMapCube.from_file_list(file_list)
    pmap_cube.plot(prefix='Original')
    # Plot the polarization map cube after the alignment.
    file_list = pipeline.file_list('stokesalign', 'pmapcube')
    pmap_cube = xBinnedPolarizationMapCube.from_file_list(file_list)
    pmap_cube.plot(prefix='Aligned')

def save_arrows_to_ds9():
    """Save the polarization arrows to a ds9 region file.
    """
    file_list = pipeline.file_list('pmapcube')
    pmap_cube = xBinnedPolarizationMapCube.from_file_list(file_list)
    output_region_file_path = os.path.join(IXPEOBSSIM_DATA, 'casa_pdpa_ds9region')
    pmap_cube.save_to_ds9(output_region_file_path, num_sigma=2., intensity_percentile=0.)

def run():
    """Run the full analysis pipeline.
    """
    simulate()
    bin_()
    display()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_casa')

# Copyright (C) 2018--2022, the ixpeobssim team.
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

import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.binning.misc import xBinnedMap
from ixpeobssim.binning.polarization import xBinnedPolarizationMapCube, xBinnedMDPMapCube


DURATION = 1000000.


def simulate():
    """Run the simulation.
    """
    pipeline.xpobssim(duration=DURATION)

def bin_():
    """Create the count map.
    """
    file_list = pipeline.file_list()
    pipeline.xpbin(*file_list, algorithm='CMAP')
    pipeline.xpbin(*file_list, algorithm='PMAPCUBE', npix=20, ebins=2)
    pipeline.xpbin(*file_list, algorithm='MDPMAPCUBE', npix=10, ebins=2)

def display():
    """
    """
    file_list = pipeline.file_list('cmap')
    map_ = xBinnedMap.from_file_list(file_list)
    pipeline.figure('count map')
    map_.plot()
    file_list = pipeline.file_list('pmapcube')
    pmap_cube = xBinnedPolarizationMapCube.from_file_list(file_list)
    pmap_cube.plot()
    file_list =  pipeline.file_list('mdpmapcube')
    mdpmap_cube = xBinnedMDPMapCube.from_file_list(file_list)
    mdpmap_cube.plot()

def run():
    """
    """
    simulate()
    bin_()
    display()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_disk')

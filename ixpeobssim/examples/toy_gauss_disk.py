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


DURATION = 100000.


def simulate():
    """Run the simulation.
    """
    pipeline.xpobssim(duration=DURATION)


def bin_():
    """Create the count map.
    """
    file_list = pipeline.file_list()
    pipeline.xpbin(*file_list, algorithm='CMAP')


def display():
    """
    """
    file_list = pipeline.file_list('cmap')
    map_ = xBinnedMap.from_file_list(file_list)
    pipeline.figure('count map')
    map_.plot()


def run():
    """
    """
    simulate()
    bin_()
    display()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('toy_gauss_disk')

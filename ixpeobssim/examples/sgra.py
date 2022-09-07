#!/usr/bin/env python
#
# Copyright (C) 2016--2018, the ixpeobssim team.
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


import os

from ixpeobssim import IXPEOBSSIM_CONFIG_REG
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
from ixpeobssim.binning.misc import xBinnedMap
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.core.spline import xInterpolatedBivariateSplineLinear
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.srcmodel.img import xFITSImage
import ixpeobssim.config.sgra as input_model
from ixpeobssim.utils.astro import read_ds9
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color

DURATION = 1000000
E_BINNING = [2., 4., 10.]
REG_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG_REG, 'sgra_analysis.reg')
RA = input_model.ROI_MODEL.ra
DEC = input_model.ROI_MODEL.dec


def generate():
    """Generate the event lists.
    """
    pipeline.xpobssim(duration=DURATION)


def make_cmap():
    """Create the Stokes cubes.
    """
    file_list = pipeline.file_list()
    pipeline.xpbin(*file_list, algorithm='CMAP', ebinalg='LIST',
                   ebinning=E_BINNING)


def plot_cmap():
    """Open files and run the standard plot.
    """
    file_list = pipeline.file_list('cmap')
    cmap = xBinnedMap.from_file_list(file_list)
    plt.figure('Counts map')
    fig = cmap.plot()
    label = 'Sgr A molecular clouds (IXPE %d Ms)' % (DURATION / 1000000.)
    #fig.add_label(0.5, 0.05, label, relative=True, size='x-large',
    #              color='white', horizontalalignment='center')
    regions = read_ds9(REG_FILE_PATH)
    ra, dec, rad = regions[0].center.ra.deg, regions[0].center.dec.deg, regions[0].radius
    #fig.show_circles(ra, dec, rad, edgecolor='green')


def select_and_bin():
    """
    """
    file_list = pipeline.file_list()
    pipeline.xpbin(*file_list, algorithm='PCUBE', ebinalg='LIST',
                   ebinning=E_BINNING)
    logger.info('Opening region file %s...' % REG_FILE_PATH)
    regions = read_ds9(REG_FILE_PATH)
    logger.info('%d region(s) found...' % len(regions))
    ra, dec, rad = regions[0].center.ra.deg, regions[0].center.dec.deg, regions[0].radius
    rad *= 60.
    logger.info('Analyzing region at ra = %s, dec = %s' % (ra, dec))
    sel_file_list = pipeline.xpselect(*file_list, ra=ra, dec=dec, rad=rad,
                                        suffix=pipeline.suffix('reg'))
    pipeline.xpbin(*sel_file_list, algorithm='PCUBE', ebinalg='LIST',
                    ebinning=E_BINNING)

def view():
    """
    """
    file_list = pipeline.file_list('reg', 'pcube')
    pcube = xBinnedPolarizationCube.from_file_list(file_list)
    pcube.plot()


def run():
    """
    """
    generate()
    make_cmap()
    plot_cmap()
    select_and_bin()
    view()


if __name__ == '__main__':
    pipeline.bootstrap_pipeline('sgra')

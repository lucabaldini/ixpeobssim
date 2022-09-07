# Copyright (C) 2022, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public Licensese as published by
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

"""Sample simulation and analysis pipeline for the SoftwareX example.

This will run a simulation for the specified configuration file, select
data for the source and the background regions, and perform a model-independent
polarization analysis in a series of energy bins.
"""
import os

import numpy

from ixpeobssim import IXPEOBSSIM_CONFIG
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.binning.polarization import xBinnedPolarizationCube

# Basic simulation parameters.
CFG_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG, 'toy_softwarex.py')
DURATION = 2000000.
# Global flag---toggle this not to overwrite existing files.
OVERWRITE = True
# Region selection: the source is a circular patch, while the background is
# a larger annulus centered in the same position (by default the reference
# position in the WCS of the original event file). All radii are in arcmin.
SRC_RAD = 0.75
BKG_INNER_RAD = 1.5
BKG_OUTER_RAD = 3.
# Energy binning for the polarization cubes.
ENERGY_BINNING = numpy.array([2., 4., 6., 8.])

# Run the simulation.
file_list = pipeline.xpobssim(configfile=CFG_FILE_PATH, duration=DURATION,
    overwrite=OVERWRITE)

# Select the source and the background regions. Note this will keep track of the
# area for each selection by setting the BACKSCAL header keyword in the output file.
src_file_list = pipeline.xpselect(*file_list, rad=SRC_RAD, suffix='src',
    overwrite=OVERWRITE)
bkg_file_list = pipeline.xpselect(*file_list, innerrad=BKG_INNER_RAD,
    rad=BKG_OUTER_RAD, suffix='bkg', overwrite=OVERWRITE)

# Create the polarization cubes.
kwargs = dict(algorithm='PCUBE', ebinalg='LIST', ebinning=ENERGY_BINNING,
    overwrite=OVERWRITE)
src_pcube_file_list = pipeline.xpbin(*src_file_list, **kwargs)
bkg_pcube_file_list = pipeline.xpbin(*bkg_file_list, **kwargs)

# Read back the polarization cubes and perform the background subtraction.
src_pcube = xBinnedPolarizationCube.from_file_list(src_pcube_file_list)
bkg_pcube = xBinnedPolarizationCube.from_file_list(bkg_pcube_file_list)
bkg_pcube *= src_pcube.backscal() / bkg_pcube.backscal()
src_pcube -= bkg_pcube

# You are good to go!
print('Polarization degree:', src_pcube.PD)
print('Polarization degree error:', src_pcube.PD_ERR)
print('Polarization angle:', src_pcube.PA, 'deg')
print('Polarization angle error:', src_pcube.PA_ERR, 'deg')

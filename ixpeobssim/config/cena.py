#!/usr/bin/env python
#
# Copyright (C) 2015--2020, the ixpeobssim team.
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

from ixpeobssim.srcmodel.roi import xChandraObservation, xChandraROIModel
from ixpeobssim.srcmodel.bkg import xPowerLawInstrumentalBkg
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.utils.astro import read_ds9
from ixpeobssim.core.fitsio import xFITSImageBase
from ixpeobssim import IXPEOBSSIM_CONFIG_FITS, IXPEOBSSIM_CONFIG_REG


__model__ = file_path_to_model_name(__file__)

# This is the input event file taken from the Chandra database (obs id 8489)
# http://cda.harvard.edu/pop/mainEntry.do
# ftp://cdaftp.cfa.harvard.edu//pub/science/ao08/cat7/8489/primary/acisf08489N003_evt2.fits.gz
EVENT_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'cena_evt.fits')

# cena_jet+core contains the definition of jet and core in the Chandra map
REG_SOURCE_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG_REG, 'cena_jet+core.reg')

regions = read_ds9(REG_SOURCE_FILE_PATH)

ROI_MODEL = xChandraROIModel(EVENT_FILE_PATH, acis='I')

pol_deg = constant(0.)
pol_ang = constant(0.)

cena = xChandraObservation('cena', pol_deg, pol_ang)
ROI_MODEL.add_source(cena)

#Add also the instrumental background
bkg = xPowerLawInstrumentalBkg()
ROI_MODEL.add_source(bkg)


def display():
    """Display the Chandra observation used to make the simulation
    and the re
    """
    IMAGE_FITS_PATH = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'cena_img.fits')
    plt.figure('Cen A, Chandra')
    image = xFITSImageBase(IMAGE_FITS_PATH)
    fig2 = image.plot(stretch='log', vmax=300)
    ax1 = plt.gca()
    ax1.annotate('Chandra 95.2 ks', xy=(0.1, 0.9), \
                 xycoords='axes fraction', fontsize = 18,\
                 color='white')



if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

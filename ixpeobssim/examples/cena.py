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

from __future__ import print_function, division


"""
To simulate Cen A we exploit the capability of the xpobssim tool to process the
provided Chandra photon list and produce an IXPE observation simulation.

The major advantage of this technique is to preserve the full correlation
between the morphology and the energy spectrum. The main simulation steps are:

* Chandra measured energies, times and positions taken as the MC truth.
* Events are down-sampled and smeared with the IXPE response functions.
* Photoelectrons angular distribution generated assuming an unpolarized source.

The input Chandra event file is taken from the Chandra database (obs id 8489)
and to obatain the spatial distribution of Cen A we used the Chandra image,
with the jet and core regions superimposed.
"""


import os

from astropy.io import fits
from astropy.wcs import WCS
import numpy

from ixpeobssim import IXPEOBSSIM_CONFIG_FITS, IXPEOBSSIM_CONFIG_REG, IXPEOBSSIM_DATA
from ixpeobssim.binning.misc import xBinnedMap
from ixpeobssim.core.fitsio import xFITSImageBase
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.evt.mdp import mdp99
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.irf import load_psf, load_modf, DEFAULT_IRF_NAME
from ixpeobssim.utils.astro import read_ds9
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt


"""Script-wide simulation and analysis settings.
"""
# This is the Chandra image of Cen A
# ftp://cdaftp.cfa.harvard.edu//pub/science/ao08/cat7/8489/primary/acisf08489N003_cntr_img2.fits.gz
IMAGE_FITS_PATH = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'cena_img.fits')

# cena_regions contains the jet, core and ULX regions in the IXPE map
REGS_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG_REG, 'cena_regions.reg')

# Source position
RA = 201.38912
DEC = -43.004776

# IXPE simulation setting
E_BINNING = [2., 8.]
DURATION = 1500000

# IXPE IRF
IRF_NAME = DEFAULT_IRF_NAME
PSF = load_psf(IRF_NAME, 1)

MODF = load_modf(IRF_NAME, 1)

REGION_LIST = ['jet', 'ulx', 'core']

def generate():
    """ Generate the simulation data starting from the Chandra event file.
    """
    pipeline.xpobssim(duration=DURATION)

def select():
    """Select the regions of cena.
    """
    file_list = pipeline.file_list()
    for region in REGION_LIST:
        region_file = os.path.join(IXPEOBSSIM_CONFIG_REG, 'cena_%s.reg' % region)
        pipeline.xpselect(*file_list, regfile=region_file, suffix=region)

def cmap():
    """Create the IXPE count map.
    """
    file_list = pipeline.file_list()
    pipeline.xpbin(*file_list, algorithm='CMAP', npix=512, pixsize=1.25)

def analyze():
    """
    Analyze the resulting simulation and calculate the MDP in the three regions
    defined in REGION_LIST.
    """
    _line = 'Results from the analysis\n'
    for region in REGION_LIST:
        energy = numpy.array([])
        src_id = numpy.array([])
        file_list = pipeline.file_list(region)
        for filepath in file_list:
            event = xEventFile(filepath)
            energy = numpy.append(energy, event.energy_data())
            src_id = numpy.append(src_id, event.srcid_data())
        _title = 'Results for %s:\n' % region
        logger.info(_title)
        _line+=_title
        ene_mask = (energy > E_BINNING[0])*(energy < E_BINNING[1])
        bkg_mask = ene_mask * (src_id == 1)  # Instrumental Background
        src_mask = ene_mask * (src_id == 0)  # All source photons
        _energy = energy[ene_mask]
        cnts_src = len(energy[src_mask])
        cnts_bkg = len(energy[bkg_mask])
        # Calculate the mdp and print the results
        mdp = mdp99(MODF.weighted_average(_energy), cnts_src, cnts_bkg)
        _fmt = '%.2f--%.2f keV: %d src counts (%.1f%%) in %d s, MDP %.2f%%\n'
        _data = (E_BINNING[0], E_BINNING[1], cnts_src,
                 100 * cnts_src / float(cnts_src + cnts_bkg), DURATION, 100 * mdp)
        _line += _fmt % _data
    logger.info(_line)


def plot(draw_regions=True):
    """ Plot the Chandra and IXPE counts map.
    """
    # Plot the IXPE map
    cmap_list = pipeline.file_list('cmap')
    full_map = xBinnedMap.from_file_list(cmap_list)
    plt.figure('Cen A, IXPE')
    fig = full_map.plot(stretch='log', zlabel='Counts')
    full_map.fits_image.recenter(RA - 0.01, DEC - 0.015, 225.)
    ax = plt.gca()
    ax.annotate('Centaurus A (IXPE 1.5 Ms)', xy=(0.1, 0.9), xycoords='axes fraction',
                fontsize=18, color='white')

    # If draw_regions flag is set to True draw the regions with labels
    if draw_regions:
        cmap_file = fits.open(cmap_list[0])
        wcs = WCS(cmap_list[0])
        regs = read_ds9(REGS_FILE_PATH)
        for i, reg in enumerate(regs):
            """SkyRegion objects currently don’t have an as_artist() or
            plot() method.
            To plot them, we need convert them to a pixel region first.
            """
            pixel_region = reg.to_pixel(wcs)
            pixel_region.plot(ax=ax)

    # Plot the original Chandra map
    plt.figure('Cen A, Chandra')
    image = xFITSImageBase(IMAGE_FITS_PATH)
    fig2 = image.plot(stretch='log', vmin=1., zlabel='Counts')
    image.recenter(RA - 0.01, DEC - 0.015, 225.)
    ax1 = plt.gca()
    ax1.annotate('Chandra 95.2 ks', xy=(0.1, 0.9), xycoords='axes fraction',
                 fontsize=18, color='white')

    if draw_regions:
        from ixpeobssim.config.cena import REG_SOURCE_FILE_PATH
        chandra_file = fits.open(IMAGE_FITS_PATH)
        chandra_wcs = WCS(IMAGE_FITS_PATH)
        r = read_ds9(REGS_FILE_PATH)
        for i, reg in enumerate(r):
            pixel_region = reg.to_pixel(chandra_wcs)
            pixel_region.plot(ax=ax1)
    plt.show()

def run():
    """
    """
    generate()
    select()
    cmap()
    analyze()
    plot(draw_regions=True)



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('cena')

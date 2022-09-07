#!/usr/bin/env python
#
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

"""xpexposure app.
"""

from __future__ import print_function, division

from astropy.io import fits
import numpy
from scipy import ndimage

from ixpeobssim.binning.exposure import xBinnedLivetimeCube, xExposureCube, create_psf_kernel
from ixpeobssim.binning.misc import xBinnedMap
from ixpeobssim.irf import load_psf, load_arf, load_vign
from ixpeobssim.irf.caldb import irf_file_path
from ixpeobssim.irf.ebounds import ENERGY_GRID, ENERGY_LO, ENERGY_HI
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger


#pylint: disable=invalid-name, too-many-locals, no-member


__description__ = \
"""Read a series of binned livetime cubes and perform the PSF convolution,
exposure calculation and correction for the effective area and/or the
modulation response functions.

The modified response files can be used, e.g., for spectro-polarimetric fitting
within XSPEC. (Note ixpeobssim offers the xpancrkey.py script to set the
ANCRFILE header keywords in the event files.)
"""

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_argument('--cmapfiles', type=str, required=True, nargs='+',
    help='path(s) to the CMAP file(s)')
PARSER.add_irfname()
PARSER.add_argument('--irftype', choices=['arf', 'mrf', 'both'], default='both',
    help='type of irf to correct and save', type=str)
PARSER.add_argument('--expunits', choices=[None, 'cm2s', 's'], default=None,
 help='flag to save the exposure cube using the required units')
PARSER.add_overwrite()


def _calculate_exposure(file_path, cmap_file, new_irfs, **kwargs):
    """Calculate the exposure and write corrected response files.
    """
    # Load the relevant objects
    ltcube = xBinnedLivetimeCube(file_path)
    cmap = xBinnedMap(cmap_file).fits_image.data
    assert ltcube.map_shape() == cmap.shape
    ontime = ltcube._img_header['ONTIME']
    du_id = int(ltcube._img_header['DETNAM'][-1])
    # Load the psf, arf and vignetting
    irf_name = kwargs.get('irfname')
    psf = load_psf(irf_name, du_id=du_id)
    arf = load_arf(irf_name, du_id=du_id)
    vign = load_vign(irf_name, du_id=du_id)
    # Create the kernel based on the psf
    kernel = create_psf_kernel(psf, ltcube.pixel_size())
    # Convolve the livetime cube with the psf
    logger.info('Applying the convolution with the PSF...')
    c_lvt = ndimage.convolve(ltcube.LIVETIME, kernel)
    theta_center = ltcube.theta_centers()
    # Create an empty array to store the exposure in several energy layers
    expo = numpy.zeros((len(ENERGY_GRID), *ltcube.map_shape()))
    # This option is used to decide whether to save or not the exposure cube
    # and which units to use (s or cm2*s).
    expunits = kwargs.get('expunits', None)
    # Loop over the energy layers and calculate the exposure. This is done by
    # multiplying the convolved livetime cube (function of theta and sky
    # coordinates) by the vignetting and summing over the theta dimension. The
    # output exposure cube will be function of energy and sky coordinates. If
    # requested the exposure is also multiplied by the effective area so that
    # the output units are cm2*s.
    for i, energy in enumerate(ENERGY_GRID):
        expo[i] = (c_lvt.T * vign(energy, theta_center)).sum(2).T
        if expunits == 'cm2s':
            expo[i] *= arf(energy)
    # Create the xExposureCube object and save it to file if requested
    expo = xExposureCube(expo, (ENERGY_LO, ENERGY_HI), expunits)
    if expunits is not None:
        logger.info('Selected units for the exposure cube: %s', expunits)
        header = ltcube.wcs.to_header()
        outfile = file_path.replace('ltcube.fits', 'expcube.fits')
        logger.info('Saving the exposure cube to %s...', outfile)
        expo.write(outfile, header, kwargs.get('overwrite'))
    # Calculate the correction factor vs energy to be applied to the effective
    # area. This is done by averaging the exposure with the CMAP, summing over
    # the sky coordinates and normalizing for the total ontime.
    aeff_corr = numpy.sum(expo.exposure * cmap / (cmap.sum() * ontime), axis=(1,2))
    # If the exposure cube units were cm2*s divide by the effective area.
    if expunits == 'cm2s':
        aeff_corr /= arf.y
    _mask = (ENERGY_GRID > 2.) * (ENERGY_GRID < 8.)
    avg_corr = numpy.mean(aeff_corr[_mask])
    logger.info('Average correction factor [2-8 keV]: %.3f', avg_corr)
    # Create the new arf and mrf reading the corresponding irf from the CALDB
    # and applying the correction factor as a function of energy.
    irf_type = [kwargs.get('irftype')]
    if irf_type[0] == 'both':
        irf_type = ['arf', 'mrf']
    for _type in irf_type:
        irf_file = irf_file_path(irf_name, du_id, _type, None, True)
        irf = fits.open(irf_file)
        irf['SPECRESP'].data.SPECRESP*= aeff_corr
        outfile = file_path.replace('.fits', '.%s' % _type)
        logger.info('Saving the new %s to %s...', _type, outfile)
        irf.writeto(outfile, overwrite=True)
        new_irfs[_type].append(outfile)
    return new_irfs


def xpexposure(**kwargs):
    """Run the app.
    """
    file_list = kwargs.get('filelist')
    cmap_files = kwargs.get('cmapfiles')
    assert len(file_list) == len(cmap_files)
    new_irfs = dict(arf=[], mrf=[])
    for file_path, cmap_file in zip(file_list, cmap_files):
        new_irfs = _calculate_exposure(file_path, cmap_file, new_irfs, **kwargs)
    return [*new_irfs['arf'], *new_irfs['mrf'], *new_irfs['mrf']]


def main():
    """main() entry point.
    """
    xpexposure(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()

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

import os

import numpy

from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.binning.polarization import xBinnedCountSpectrum
from ixpeobssim.core.spline import xUnivariateSpline
<<<<<<< HEAD
from ixpeobssim.instrument.gpd import FIDUCIAL_AREA
from ixpeobssim.instrument.mma import FIDUCIAL_BACKSCAL
=======
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.instrument.gpd import GPD_PHYSICAL_AREA
from ixpeobssim.instrument.mma import FOCAL_LENGTH
>>>>>>> main
from ixpeobssim.irf.ebounds import channel_to_energy, ENERGY_STEP
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, residual_plot

'''
Create a background template model starting from a series of PHA1
background files.

The PHA1 files should be prepared from a dark field with a arbitrary
shapes that have been cut from one or more files and have a BACKSCAL keyword
defined.

This is suitable both for residual and total background, depending on the user
needs.

The count spectrum, normalized by the backscal, the fiducial area of the
detector and by the bin width, is parametrized with a non interpolated
spline and written to file to be used later with an appropriate config file.

The output file is written on a regular energy grid as a simple text file
with two columns---energy and background rate.
'''


def smooth_PDF(PDF, artificial_grad = 1.e-6):
    ''' Replaces the initial range (E<1keV) of the spline with a monothonic
    function and adjusts the interval of negative values with an artificial
    gradient.

    Note that for very dramatically compromised spectral shapes this does not,
    nor is intended to work, and we should work on the spline smoothing 
    parameter instead. 
    '''
    limit = 17 # E< keV
    PDF[0:limit] = numpy.sort(PDF[0:limit])
    is_zero_idx = numpy.where (PDF[0:limit]<=0)
    max_zero = numpy.max(is_zero_idx)
    # turn negative numbers into zeros
    for j in range (limit):
        if PDF[j]<0: PDF[j]=0
    #We don't want to approach zero with a large derivative
    for j in reversed(range(max_zero+1)):
        PDF[j] = PDF[j+1]-artificial_grad
    return(PDF)


def create_backgound_template(phalist, ssmooth, outfile, emin=0.01):
    '''The core function: takes a list of pha1 files and creates the background \
        template. The parameters are taken from input.
    '''


    logger.info (f'loading background spectra from {phalist}...')

    # Loop over all input background files:
    # Load the raw count spectrum and convert PI channels in keV.
    # Calculate the scaling factors and convert in proper units
    livetime_total = 0
    for file in phalist:
        spec = xBinnedCountSpectrum.from_file_list([file])
        # Load the livetime and divide all quantities for weighted average
        livetime = spec.spectrum_header['LIVETIME']
        backscal = spec.spectrum_header['BACKSCAL']
        # Rate * livetime = total signal
        spec.RATE *= livetime
        spec.STAT_ERR *= livetime
        livetime_total += livetime
        # Divide by the fraction of the total area
        logger.info('Correcting for the extraction radius...')
        area_frac = backscal / FIDUCIAL_BACKSCAL
        logger.info (f'Fraction of the total area: {area_frac}')
        spec.RATE /= area_frac
        spec.STAT_ERR /= area_frac
        if 'avg_spec' in locals():
            avg_spec += spec
        else:
            avg_spec = spec
        plt.figure('Single file input spectra vs average')
        plt.semilogy(channel_to_energy(spec.CHANNEL),spec.RATE / livetime, lw = 1)
    plt.semilogy(channel_to_energy(avg_spec.CHANNEL),avg_spec.RATE / livetime_total,
        linewidth = 4, label = 'mean')
    plt.xlabel('Energy [keV]')
    plt.ylabel('Background spectrum')
    plt.grid(which='both')
    plt.legend()

<<<<<<< HEAD
    # Scale the single object for the overall livetime
    # Integrated signal / total livetime = rate again
    logger.info (f'scaling back for the total livetime of {livetime_total}')
    avg_spec.RATE /= livetime_total
    avg_spec.STAT_ERR /= livetime_total
    # To physical units, also accounted for the backscal division
    logger.info('Converting into physical units...')
    scale = 1. / (FIDUCIAL_AREA / 100.) / ENERGY_STEP
    logger.info('Region backscal: %.3e', scale)
    flux = avg_spec.RATE * scale
    flux_err = avg_spec.STAT_ERR * scale
    energy = channel_to_energy(avg_spec.CHANNEL)
=======
    emin : float
        The minimum energy for the spline.
    """
    # The file name, at this point, is hard-coded---this can be made more
    # user-friendly in the future.
    file_name = 'ixpe01903701_det%d_evt2_v02_bkg_pha1.fits'
    file_list = [os.path.join(IXPEOBSSIM_BKG_DATA, file_name % du_id) for du_id in DU_IDS]
    output_file_name = 'bkg_smcx1_01903701.txt'
    output_file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'ascii', output_file_name)

    # Load the raw count spectrum and convert PI channels in keV. Note we are
    # divding by the number of detector units.
    spec = xBinnedCountSpectrum.from_file_list(file_list)
    num_detectors = len(file_list)
    energy = channel_to_energy(spec.CHANNEL)
    rate = spec.RATE / num_detectors
    rate_err = spec.STAT_ERR / num_detectors

    # Calculate the scaling factor to compensate for the source cut in the
    # underlying PHA1 file, i.e., convert the extraction radius from arcmin to mm
    # and scale to the full fiducial area of the GPD.
    logger.info('Correcting for the extraction radius...')
    radius = FOCAL_LENGTH * numpy.radians(arcmin_to_degrees(extraction_radius))
    scale = 1. / (1. - numpy.pi * radius**2. / GPD_PHYSICAL_AREA)
    logger.info('Scaling factor: %.3f', scale)
    rate *= scale
    rate_err *= scale

    # Convert in the proper units, i.e., from s^{-1} to cm^{-1} s^{-1} keV^{-1}.
    # This is done by divding by the detector fiducial area (converted from mm^2
    # to cm^2) and the 40 eV step of the energy grid used to create the response
    # files.
    logger.info('Converting into physical units...')
    scale = 1. / (GPD_PHYSICAL_AREA / 100.) / ENERGY_STEP
    logger.info('Scaling factor: %.3e', scale)
    flux = rate * scale
    flux_err = rate_err * scale

>>>>>>> main
    # Create a non-interpolating spline---note that we are cutting at a minimum
    # energy to avoid the peak in channel 0.
    mask = energy > emin
    spline = xUnivariateSpline(energy[mask], flux[mask], s=ssmooth, k=3)
    # Create the output text file with the spline data.
    logger.info('Writing output file to %s...', outfile)
    x = numpy.linspace(emin, energy.max(), 250)
    y = spline(x).clip(0.)
    y = smooth_PDF(y)
    
    with open(outfile, 'w') as output_file:
        for _x, _y in zip(x, y):
            output_file.write('%.5e   %.5e\n' % (_x, _y))
    logger.info('Done.')

    # Plotting stuff.
    ax1, ax2 = residual_plot('background template')
    plt.errorbar(energy, flux, flux_err, fmt='o', label='Background data (all DUs)')
    spline.plot(zorder=1, label='Spline approximation')
    setup_gca(ylabel='Background rate [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]', logy=True,
        grids=True, xmax=energy.max(), legend=True)
    plt.sca(ax2)
    plt.errorbar(energy, flux - spline(energy), flux_err, fmt='o')
    dy = 5.e-3
    setup_gca(xlabel='Energy [keV]', ymin=-dy, ymax=dy, grids=True, xmax=energy.max())
    plt.show()


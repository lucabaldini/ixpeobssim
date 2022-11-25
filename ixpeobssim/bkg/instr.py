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
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.instrument.gpd import GPD_PHYSICAL_AREA, GPD_DEFAULT_FIDUCIAL_HALF_SIDE_X,\
    GPD_DEFAULT_FIDUCIAL_HALF_SIDE_Y, fiducial_area
from ixpeobssim.instrument.mma import FOCAL_LENGTH, fiducial_backscal
from ixpeobssim.irf.ebounds import channel_to_energy, ENERGY_STEP
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, residual_plot



def smooth_PDF(PDF, artificial_grad = 1.e-6):
    """Replaces the initial range (E<1keV) of the spline with a monothonic
    function and adjusts the interval of negative values with an artificial
    gradient.

    Note that for very dramatically compromised spectral shapes this does not,
    nor is intended to work, and we should work on the spline smoothing
    parameter instead.
    """
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
    """Create a background template model starting from a series of PHA1 background files.

    The PHA1 files should be prepared from a dark field with a arbitrary shape that
    have been cut from one or more files and have a BACKSCAL keyword defined.

    This is suitable both for residual and total background, depending on the user
    needs.

    The count spectrum, normalized by the backscal, the fiducial area of the detector
    and by the bin width, is parametrized with a non interpolated spline and written
    to file to be used later with an appropriate config file.

    The output file is written on a regular energy grid as a simple text file with
    two columns---energy and background rate.
    """
    logger.info (f'loading background spectra from {phalist}...')

    # Loop over all input background files:
    # Load the raw count spectrum and convert PI channels in keV.
    # Calculate the scaling factors and convert in proper units
    livetime_total = 0
    half_side_x = GPD_DEFAULT_FIDUCIAL_HALF_SIDE_X
    half_side_y = GPD_DEFAULT_FIDUCIAL_HALF_SIDE_Y
    det_backscal = fiducial_backscal(half_side_x, half_side_y)
    det_area = fiducial_area(half_side_x, half_side_y)
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
        area_frac = backscal / det_backscal
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

    # Scale the single object for the overall livetime
    # Integrated signal / total livetime = rate again
    logger.info (f'scaling back for the total livetime of {livetime_total}')
    avg_spec.RATE /= livetime_total
    avg_spec.STAT_ERR /= livetime_total
    # To physical units, also accounted for the backscal division
    logger.info('Converting into physical units...')
    scale = 1. / (det_area / 100.) / ENERGY_STEP
    logger.info('Region backscal: %.3e', scale)
    flux = avg_spec.RATE * scale
    flux_err = avg_spec.STAT_ERR * scale
    energy = channel_to_energy(avg_spec.CHANNEL)
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

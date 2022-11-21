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

"""Facilities to create templates for the instrumental background.
"""

from __future__ import print_function, division
from operator import is_

import os

import numpy

from ixpeobssim import IXPEOBSSIM_BKG_DATA, IXPEOBSSIM_SRCMODEL
from ixpeobssim.binning.polarization import xBinnedCountSpectrum
from ixpeobssim.core.spline import xUnivariateSpline
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.instrument.gpd import FIDUCIAL_AREA
from ixpeobssim.instrument.mma import FOCAL_LENGTH
from ixpeobssim.instrument.mma import FIDUCIAL_BACKSCAL
from ixpeobssim.irf.ebounds import channel_to_energy, ENERGY_STEP
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, residual_plot
from ixpeobssim.utils.units_ import arcmin_to_degrees
from glob import glob
from astropy.io import fits



def create_backgound_template(spline_smoothing=5.e-4, emin=0.01):
    """Create a background template model starting from a series of PHA1
    background files.

    The PHA1 files are intented to be prepared starting from a field with a
    point sources, that has been preliminary been removed via a circular cut
    with a specified radius---the latter is used internally to reascale the
    rate to the full detector surface.

    The count spectrum, normalized by the backscal, the fiducial area of the 
    detector and by the bin width, is then parametrized with a non interpolated
    spline and written to file to be used later.

    The output file is written on a regular energy grid as a simple text file
    with two columns---energy and background rate.

    Args
    ----
    spline_smoothing : float
        The smoothing coefficient ("s" argument in the scipy documentation) used
        for the non interpolating spline. Note this is very important, as it
        controls the level at which the spline is capturing the fluctuations of
        the input data points. (s=0 is effectively an interpolating spline, but
        the actual value depends on the scale of the input data and it's not
        trivial to establish a priori.)

    emin : float
        The minimum energy for the spline.
    """

    # The file name, at this point, is hard-coded---this can be made more
    # user-friendly in the future.
    #file_name = 'ixpe01005401_det%d_evt2_v02_clean_bkg_pha1.fits'
    #file_list = [os.path.join(IXPEOBSSIM_BKG_DATA, file_name % du_id) for du_id in DU_IDS]
    file_list = glob(f'{IXPEOBSSIM_BKG_DATA}/*5701*pha1.fits')
    logger.info (f'loading background spectra from {file_list}...')
    output_file_name = 'bkg_conservative_3c279.txt'
    output_file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'ascii', output_file_name)

    # Loop over all input background files:
    # Load the raw count spectrum and convert PI channels in keV. 
    # Calculate the scaling factors and convert in proper units
    livetime_total = 0
    for file in file_list:
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
        plt.loglog(channel_to_energy(spec.CHANNEL),spec.RATE / livetime, linewidth = 1)
    plt.loglog(channel_to_energy(avg_spec.CHANNEL),avg_spec.RATE / livetime_total, linewidth = 4, label = f'mean')
    plt.xlabel('Energy [keV]')
    plt.ylabel('Background spectrum')
    plt.legend()
    plt.show()

    
    # Scale the single object for the overall livetime
    # Integrated signal / total livetime = rate again
    logger.info (f'scaling back for the total livetime of {livetime_total}')
    avg_spec.RATE /= livetime_total
    avg_spec.STAT_ERR /= livetime_total
    # To physical units, also accounted for the backscal division
    logger.info('Converting into physical units...')
    scale = 1. / (FIDUCIAL_AREA / 100.) / ENERGY_STEP
    logger.info('Region backscal: %.3e', scale)
    rate = avg_spec.RATE
    rate_err = avg_spec.STAT_ERR
    flux = rate * scale
    flux_err = rate_err * scale
    energy = channel_to_energy(avg_spec.CHANNEL)
    



    #spec = xBinnedCountSpectrum.from_file_list(file_list)
    #livetime = spec.spectrum_header['LIVETIME']
    #num_detectors = len(file_list)
    #energy = channel_to_energy(spec.CHANNEL)
    #rate = spec.RATE / num_detectors
    #rate_err = spec.STAT_ERR / num_detectors

    # Calculate the scaling factor to compensate for the source cut in the
    # underlying PHA1 file, i.e., convert the extraction radius from arcmin to mm
    # and scale to the full fiducial area of the GPD.
    #logger.info('Correcting for the extraction radius...')
    #radius = FOCAL_LENGTH * numpy.radians(arcmin_to_degrees(extraction_radius))
    #scale = 1. / (1. - numpy.pi * radius**2. / FIDUCIAL_AREA)
    #logger.info('Scaling factor: %.3f', scale)
    #rate *= scale
    #rate_err *= scale

    # Convert in the proper units, i.e., from s^{-1} to cm^{-1} s^{-1} keV^{-1}.
    # This is done by divding by the detector fiducial area (converted from mm^2
    # to cm^2) and the 40 eV step of the energy grid used to create the response
    # files.
    #logger.info('Converting into physical units...')
    #scale = 1. / (FIDUCIAL_AREA / 100.) / ENERGY_STEP
    #logger.info('Scaling factor: %.3e', scale)
    #flux = rate * scale
    #flux_err = rate_err * scale


    # Create a non-interpolating spline---note that we are cutting at a minimum
    # energy to avoid the peak in channel 0.
    mask = energy > emin
    spline = xUnivariateSpline(energy[mask], flux[mask], s=spline_smoothing, k=3)
    # Create the output text file with the spline data.
    logger.info('Writing output file to %s...', output_file_path)
    x = numpy.linspace(emin, energy.max(), 250)
    y = spline(x).clip(0.)
    with open(output_file_path, 'w') as output_file:
        for _x, _y in zip(x, y):
            output_file.write('%.5e   %.5e\n' % (_x, _y))
    logger.info('Done.')

    # Plotting stuff.
    ax1, ax2 = residual_plot('background template')
    plt.errorbar(energy, flux, flux_err, fmt='o', label='Background data (all DUs)')
    spline.plot(zorder=3, label='Spline approximation')
    setup_gca(ylabel='Background rate [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]', logy=True,
        grids=True, xmax=energy.max(), legend=True)
    plt.sca(ax2)
    plt.errorbar(energy, flux - spline(energy), flux_err, fmt='o')
    dy = 5.e-3
    setup_gca(xlabel='Energy [keV]', ymin=-dy, ymax=dy, grids=True, xmax=energy.max())
    plt.show()



if __name__ == '__main__':
    create_backgound_template()

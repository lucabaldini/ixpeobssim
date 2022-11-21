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

from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.binning.polarization import xBinnedCountSpectrum
from ixpeobssim.core.spline import xUnivariateSpline
from ixpeobssim.instrument.gpd import FIDUCIAL_AREA
from ixpeobssim.instrument.mma import FIDUCIAL_BACKSCAL
from ixpeobssim.irf.ebounds import channel_to_energy, ENERGY_STEP
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, residual_plot


__description__ = \
'''
Create a background template model starting from a series of PHA1
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
'''

parser = xArgumentParser(description=__description__)
parser.add_argument('phalist', nargs='+', help='path(s) to the pha1 \
                    spectrum file(s) upon which to build the template')
parser.add_argument('--ssmooth', type=float, default=5.e-4,
    help='The smoothing coefficient ("s" argument in the scipy documentation) \
        used for the non interpolating spline. Note this is very important, as \
        it controls the level at which the spline is capturing the fluctuations \
        of the input data points. (s=0 is effectively an interpolating spline, \
        but the actual value depends on the scale of the input data and it is \
        not trivial to establish a priori.)')
parser.add_outfile(default=os.path.join(IXPEOBSSIM_SRCMODEL, 'ascii',
        'instrumental_bkg_template.txt'))

def create_backgound_template(emin=0.01, **kwargs):
    '''The core function: takes a list of pha1 files and creates the background \
        template. The parameters are taken from input.
    '''

    phalist = kwargs.get('phalist')
    ssmoth = kwargs.get('ssmooth')
    outfile = kwargs.get('outfile')
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
    # Create a non-interpolating spline---note that we are cutting at a minimum
    # energy to avoid the peak in channel 0.
    mask = energy > emin
    spline = xUnivariateSpline(energy[mask], flux[mask], s=ssmoth, k=3)
    # Create the output text file with the spline data.
    logger.info('Writing output file to %s...', outfile)
    x = numpy.linspace(emin, energy.max(), 250)
    y = spline(x).clip(0.)
    with open(outfile, 'w') as output_file:
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
    args = parser.parse_args()
    create_backgound_template(**vars(args))

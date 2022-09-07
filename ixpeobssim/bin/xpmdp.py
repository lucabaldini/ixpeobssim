#!/usr/bin/env python
#
# Copyright (C) 2016--2020, the ixpeobssim team.
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


__description__ = \
"""Calculate the minimum detectable polarization (MDP) for a given source model.

The program takes the same python configuration modules that can be fed
into xpobbsim and calculates the MDP at the 99% CL by direct numerical
integration of the input spectrum (i.e., there are no event lists involved).
As a consequence, part of the richness of the detector response (most
notably, the energy dispersion and the vignetting) is not captured here.

For use cases beyond simple stationary point sources, the use of xpobbsim and
xpbin mode are recommended, as that approach offers the maximum flexibility.
"""

import os

from ixpeobssim.instrument import DU_IDS
from ixpeobssim.binning.base import xEventBinningBase
from ixpeobssim.irf import load_arf, load_modf
from ixpeobssim.srcmodel.spectrum import xCountSpectrum
from ixpeobssim.srcmodel.roi import xPeriodicPointSource
from ixpeobssim.evt.mdp import xMDPTable
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.srcmodel import import_roi


EBIN_ALGS = ['FILE', 'LIN', 'LOG', 'LIST']



def _build_source(**kwargs):
    """Build the source component for the sensitivity instrument.
    """
    roi_model = import_roi(kwargs.get('configfile'))
    return roi_model.source_by_id(kwargs.get('srcid'))


def _build_spectrum(source, du_id, **kwargs):
    """Build the count spectrum for the source and a specific DU.
    """
    aeff = load_arf(kwargs.get('irfname'), du_id)
    if isinstance(source, xPeriodicPointSource):
        grid = numpy.linspace(kwargs.get('phasemin'), kwargs.get('phasemax'), 100)
        scale_factor = kwargs.get('duration')
    else:
        grid = source.sampling_time_grid(0., kwargs.get('duration'))
        scale_factor = 1.
    return xCountSpectrum(source.photon_spectrum, aeff, grid, source.column_density,
                          source.redshift, scale_factor=scale_factor)


def calculate_mpd(source_factory, spectrum_factory, **kwargs):
    """Actual routine calculating the MDP.

    This is setup in this fashion so that the same code can be reused by
    xppimms, that is doing a very similar calculation, with minor differences as
    to how the underlying source and spectra are calculated.
    """
    # Build the source.
    source = source_factory(**kwargs)
    logger.info(source)
    # Cache a few things.
    duration = kwargs.get('duration')
    deadtime = kwargs.get('deadtime')
    eef = kwargs.get('eef')
    ebinning = xEventBinningBase.make_energy_binning(**kwargs)
    # Create an empty MDP table.
    mdp_table = xMDPTable.empty(duration, ebinning)
    # Start the loop over the three detector units.
    for du_id in DU_IDS:
        modf = load_modf(kwargs.get('irfname'), du_id)
        count_spectrum = spectrum_factory(source, du_id, **kwargs)
        _table = count_spectrum.build_mdp_table(ebinning, modf)
        if deadtime > 0.:
            # Deadtime correction---note we are using the broadband count spectrum
            # integral (defaulting to 1--12 keV) for the deadtime estimation,
            # independently from the actual MDP table.
            correction = 1. - count_spectrum.num_expected_counts() * deadtime / duration
            logger.info('Correcting the MDP table for deadtime (correction = %.5f)', correction)
            _table.scale(correction)
        mdp_table += _table
    # Finalize the table.
    mdp_table.set_source(source)
    if eef < 1.:
        logger.info('Correcting the MDP table for the PSF EEF (correction = %.5f)', eef)
        mdp_table.scale(eef)
    print(mdp_table)
    file_path = kwargs.get('outfile')
    if file_path is not None:
        mdp_table.save_ascii(file_path, **kwargs)
    return mdp_table


def xpmdp(**kwargs):
    """Calculate the MDP.
    """
    return calculate_mpd(_build_source, _build_spectrum, **kwargs)


# Command-line switches.
PARSER = xArgumentParser(description=__description__)
PARSER.add_outfile()
PARSER.add_configfile(required=True)
PARSER.add_srcid()
PARSER.add_irfname()
PARSER.add_duration(10000.)
PARSER.add_deadtime()
PARSER.add_eef()
PARSER.add_phasebounds(0., 1.)
PARSER.add_ebinning(default_emin=2., default_emax=8., ebin_algs=EBIN_ALGS,
                    default_ebinalg='LOG', default_ebins=4)
PARSER.add_overwrite()


def main():
    xpmdp(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

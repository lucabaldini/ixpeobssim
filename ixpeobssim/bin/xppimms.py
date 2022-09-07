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

from __future__ import print_function, division, absolute_import


__description__ = \
"""Calculate the minimum detectable polarization (MDP) for simple stationary
point sources with a power-law spectrum, according to the source parameters
provided through the command-line switches.

(Note the default power-law normalization and index and the default duration
correspond to the setup for the IXPE level-1 polarization requirement).

The program calculates the MDP at the 99% CL by direct numerical integration of
the input spectrum and  response functions (i.e., there are no event lists
involved). As a consequence, part of the richness of the detector response
(most notably, the energy dispersion and the effective area vignetting) is not
captured.

For use cases beyond simple stationary point sources, the use of xpobbsim,
xpselect and xpbin are recommended, as that approach offers the maximum
flexibility.
"""

from ixpeobssim.irf import load_arf
from ixpeobssim.srcmodel.spectrum import xCountSpectrum
from ixpeobssim.srcmodel.roi import xPointSource
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.bin.xpmdp import calculate_mpd, EBIN_ALGS
from ixpeobssim.utils.argparse_ import xArgumentParser


def _build_source(**kwargs):
    """Build the source for the MDP calculation.

    This returns essentially a fictional source at ra, dec = 0., 0., with
    zero polarization degree and angle, and the appropriate power-law
    spectrum.
    """
    photon_spectrum = power_law(kwargs.get('norm'), kwargs.get('index'))
    return xPointSource('Point source', 0., 0., photon_spectrum, constant(0.),
                        constant(0.), kwargs.get('nH'),
                        redshift=kwargs.get('redshift'))


def _build_spectrum(source, du_id, **kwargs):
    """Build the count spectrum for the source and a specific DU.

    Note this is simpler than the function in xpmdp, as we are not interested
    in periodic sources.
    """
    aeff = load_arf(kwargs.get('irfname'), du_id)
    grid = source.sampling_time_grid(0., kwargs.get('duration'))
    return xCountSpectrum(source.photon_spectrum, aeff, grid, source.column_density,
                          source.redshift)


def xppimms(**kwargs):
    """Calculate the MDP.
    """
    return calculate_mpd(_build_source, _build_spectrum, **kwargs)



# Command-line switches.
PARSER = xArgumentParser(description=__description__)
PARSER.add_outfile()
PARSER.add_irfname()
PARSER.add_pl_norm()
PARSER.add_pl_index()
PARSER.add_column_density()
PARSER.add_redshift()
PARSER.add_duration(864000.)
PARSER.add_deadtime()
PARSER.add_eef()
PARSER.add_ebinning(default_emin=2., default_emax=8., ebin_algs=EBIN_ALGS,
                    default_ebinalg='LOG', default_ebins=4)
PARSER.add_overwrite()



def main():
    """Main function.
    """
    args = PARSER.parse_args()
    xppimms(**args.__dict__)


if __name__ == '__main__':
    main()

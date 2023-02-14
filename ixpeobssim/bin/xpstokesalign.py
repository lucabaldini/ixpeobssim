#!/usr/bin/env python
#
# Copyright (C) 2019 the ixpeobssim team.
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
"""Align the reconstructed Stokes parameters to a given polarization model.

This is an application processing a photon list and `rotating` the Stokes parameters
so that, on an event-by-event basis, the zero for the measurement of the
photoelectron direction is aligned to a given input model at the position of the
event.

In the simplest form the alignment can be either radial or tangential, which is
achieved by selecting the RAD or TAN modes, respectively, and optionally passing
the right ascension and the declination of the center for the rotation via the
--ra and --dec command-line switches.

Additionally, the user can feed into tha application pair of FITS images
(in either the Q/U, X/Y components of the polarization vector or polarization
degree/angle) in the same exact fashion of the machinery used for simulating
complex polarization patterns for extended sources in ixpeobssim. This is
achieved via a combination of the QU, XY of PDA modes, along with the proper
model files passed to the --modelfiles command-line switch.

Supported alignment algorithms:
- RAD: radial
- TAN: tangential
- QU : generic polarization model (from FITS maps of U and Q)
- XY : generic polarization model (from FITS maps of polarization components)
- PDA: generic polarization angle model (from a single FITS map)

Known limitations:
- the applications is currently not supporting Stokes sky-cubes, i.e., models
  where the polarization angle is energy-dependent (although this can be
  implemented if there is need for);
- when running in PDA mode, although in principle the map of the polarization
  angle would be enough to operate the rotation, for techincal reasons having
  to do with the internals of the code, we still need to pass the maps of both
  the polarization degree and angle (and, again, this is a nuisance that we
  can overcome if there is really need for).
"""

import os

import numpy

from ixpeobssim.evt.align import align_stokes_parameters
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.evt.kislat2015 import xStokesAnalysis
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.srcmodel.polarization import xStokesSkyMap, xRadialPolarizationField,\
    xTangentialPolarizationField
from ixpeobssim.utils.argparse_ import xArgumentParser, argparse

# pylint: disable=invalid-name


ALIGN_MODES = ('RAD', 'TAN', 'QU', 'PDA')

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_argument('--mode', choices=ALIGN_MODES, type=str, required=True,
                    default=argparse.SUPPRESS,
                    help='the alignment mode')
PARSER.add_argument('--modelfiles', nargs='+',
                    help='path to the input model file(s)')
PARSER.add_argument('--ra', type=float, default=None,
                    help='RA of the model center in decimal degrees')
PARSER.add_argument('--dec', type=float, default=None,
                    help='Dec of the model center in decimal degrees')
PARSER.add_suffix(None)
PARSER.add_mc()
PARSER.add_overwrite()



def xpstokesalign(**kwargs):
    """Run the phi alignmnent.
    """
    mode = kwargs.get('mode')
    logger.info('Alignment mode: "%s"', mode)
    assert mode in ALIGN_MODES
    file_list = kwargs.get('filelist')
    outlist = []
    for file_path in file_list:
        assert file_path.endswith('.fits')
        outfile = file_path.replace('.fits', '_stokesalign.fits')
        suffix = kwargs.get('suffix')
        if suffix is not None:
            outfile = outfile.replace('.fits', '_%s.fits' % suffix)
        outlist.append(outfile)
        if os.path.exists(outfile) and not kwargs.get('overwrite'):
            logger.info('Output file %s already exists.', outfile)
            logger.info('Remove it or set "overwrite = True" to overwite it.')
            continue
        # Open the input file and recover the relevant columns.
        input_file = xEventFile(file_path)
        ra, dec = input_file.sky_position_data(kwargs.get('mc'))
        q = input_file.q_data()
        u = input_file.u_data()
        # Retrieve the position of the center for RAD/TAN alignmnent.
        obj_ra, obj_dec = input_file.wcs_reference()
        ra0 = kwargs.get('ra') or obj_ra
        dec0 = kwargs.get('dec') or obj_dec
        # Calculate the reference angle for the alignment.
        if mode == 'RAD':
            model = xRadialPolarizationField(ra0, dec0)
        elif mode == 'TAN':
            model = xTangentialPolarizationField(ra0, dec0)
        elif mode == 'QU':
            model = xStokesSkyMap.load_from_qu(*kwargs.get('modelfiles'))
        elif mode == 'PDA':
            model = xStokesSkyMap.load_from_pda(*kwargs.get('modelfiles'))
        phi0 = model.polarization_angle(ra, dec)
        # Recalculate the Stokes parameters.
        q0 = xStokesAnalysis.stokes_q(phi0, weights=None)
        u0 = xStokesAnalysis.stokes_u(phi0, weights=None)
        q, u = align_stokes_parameters(q, u, q0, u0)
        # Overwrite the proper columns in the event list.
        data = input_file.hdu_list['EVENTS'].data
        data['Q'] = q
        data['U'] = u
        # And we're ready to write the output file!
        logger.info('Writing output file %s', outfile)
        input_file.hdu_list.writeto(outfile, overwrite=True)
    return outlist


def main():
    """main() entry point
    """
    xpstokesalign(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

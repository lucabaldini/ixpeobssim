#!/usr/bin/env python
#
# Copyright (C) 2021 the ixpeobssim team.
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
"""Write a text file to be fed as an input custom spectrum to ixpesim to
simulate a realistic astrophysical spectrum.

At the moment the spectrum is limited to a power law with an adjustable index,
and the H column density can be controlled via command line. The effect of the
mirror effective area and the transparency of the UV filter are taken into
account.

The output is a simple text file with two columns---the energy in MeV and the
flux in arbitrary units.
"""

import os

import numpy

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.irfgen.du import uv_filter_transparency
from ixpeobssim.irfgen.gpd import window_transparency
from ixpeobssim.irfgen.mma import effective_area
from ixpeobssim.srcmodel.gabs import xInterstellarAbsorptionModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger


PARSER = xArgumentParser(description=__description__)
PARSER.add_ebounds(1., 12.)
PARSER.add_argument('--ebins', type=int, default=250, help='number of energy bins')
PARSER.add_pl_index()
PARSER.add_column_density()


def calculate_spectrum(energy, pl_index, column_density, du_id=1):
    """Calculate the actual spectrum, e.g., convolve the power law with the
    MMA effective area and the UV filter transparency, along with the
    interstellar absoprption if needed.

    .. warning::

       Since ixpesim is not equipped to include the long list of contaminants
       in the detector simulation (and it makes little sense to do so) we need
       to multiply the input spectrum by an energy-dependent correction factor
       of the order of unity to take into account that.
    """
    spectrum = power_law(1., pl_index)(energy)
    if column_density > 0:
        ism = xInterstellarAbsorptionModel()
        trans = ism.transmission_factor(column_density)
        spectrum *= trans(energy)
    spectrum *= effective_area(energy, du_id)
    spectrum *= uv_filter_transparency(energy)
    # Correct for the different modeling of the Be window assemblies between
    # ixpeobssim and ixpesim.
    scale = window_transparency(energy) / window_transparency(energy, None)
    spectrum *= scale
    # Normalize to the maximum value, although this is irrelevant.
    spectrum /= spectrum.max()
    return spectrum

def xpsimspec(**kwargs):
    """Dump the glorious ixpesim custom spectrum.
    """
    energy = numpy.linspace(kwargs['emin'], kwargs['emax'], kwargs['ebins'])
    pl_index = kwargs['index']
    column_density = kwargs['nH']
    spectrum = calculate_spectrum(energy, pl_index, column_density)
    file_name = 'ixpesim_spec_gamma%.2f_nh%.2e.dat' % (pl_index, column_density)
    file_path = os.path.join(IXPEOBSSIM_DATA, file_name)
    logger.info('Writing custom spectrum to %s...', file_path)
    with open(file_path, 'w') as output_file:
        for _energy, _spec in zip(energy, spectrum):
            output_file.write('%.6f %.6f\n' % (0.001 * _energy, _spec))
    logger.info('Done.')



def main():
    xpsimspec(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()

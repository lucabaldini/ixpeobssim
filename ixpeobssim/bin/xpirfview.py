#!/usr/bin/env python
#
# Copyright (C) 2015, the ixpeobssim team.
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
"""Generic display application for viewing response files.

This application allows to display the content of IXPE effective area,
point-spread-function, energy dispersion and modulation factor response files.
"""


from ixpeobssim.utils.logging_ import abort
from ixpeobssim.irf import peek_irf_type
from ixpeobssim.irf.arf import xEffectiveArea
from ixpeobssim.irf.modf import xModulationFactor
from ixpeobssim.irf.mrf import xModulationResponse
from ixpeobssim.irf.psf import xPointSpreadFunction
from ixpeobssim.irf.rmf import xEnergyDispersion
from ixpeobssim.irf.vign import xVignetting
from ixpeobssim.utils.matplotlib_ import plt, save_gcf


CLASS_DICT = {
    'arf' : xEffectiveArea,
    'modf': xModulationFactor,
    'mrf' : xModulationResponse,
    'rmf' : xEnergyDispersion,
    'psf' : xPointSpreadFunction,
    'vign': xVignetting
}


def xpirfview(file_path, save, output_folder):
    """Quick FITS image viewer.
    """
    _type = peek_irf_type(file_path)
    if not _type in CLASS_DICT.keys():
        abort('Unrecognized irf type (%s)' % _type)
    irf = CLASS_DICT[_type](file_path)
    irf.plot()
    if save:
        save_gcf(output_folder)
    plt.show()


def main():
    from ixpeobssim.utils.argparse_ import xArgumentParser
    parser = xArgumentParser(description=__description__)
    parser.add_file()
    parser.add_save()
    parser.add_outfolder()
    args = parser.parse_args()
    xpirfview(args.file, args.save, args.outfolder)



if __name__ == '__main__':
    main()

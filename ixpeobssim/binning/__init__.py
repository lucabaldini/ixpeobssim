# Copyright (C) 2015--2022, the ixpeobssim team.
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

"""Mapping of the binning data types.
"""

from __future__ import print_function, division

import ixpeobssim.binning.detector
import ixpeobssim.binning.exposure
import ixpeobssim.binning.misc
import ixpeobssim.binning.polarization


# Basic dictionary for binned product I/O---this is used, e.g., in xpbin.
BINNING_WRITE_DICT = {
    'PHA1'      : polarization.xEventBinningPHA1,
    'PHA1Q'     : polarization.xEventBinningPHA1Q,
    'PHA1U'     : polarization.xEventBinningPHA1U,
    'PHA1QN'    : polarization.xEventBinningPHA1QN,
    'PHA1UN'    : polarization.xEventBinningPHA1UN,
    'CMAP'      : misc.xEventBinningCMAP,
    'MDPMAP'    : polarization.xEventBinningMDPMAP,
    'MDPMAPCUBE': polarization.xEventBinningMDPMAPCUBE,
    'PCUBE'     : polarization.xEventBinningPCUBE,
    'PMAP'      : polarization.xEventBinningPMAP,
    'PMAPCUBE'  : polarization.xEventBinningPMAPCUBE,
    'PP'        : misc.xEventBinningPP,
    'ARMAP'     : detector.xEventBinningARMAP,
    'EFLUX'     : detector.xEventBinningEFLUX,
    'LC'        : misc.xEventBinningLC,
    'LTCUBE'    : exposure.xEventBinningLTCUBE
}


# Basic dictionary for binned product I/O---this is used, e.g., in xpbinview.
BINNING_READ_DICT = {
    'PHA1'      : polarization.xBinnedCountSpectrum,
    'PHA1Q'     : polarization.xBinnedCountSpectrum,
    'PHA1U'     : polarization.xBinnedCountSpectrum,
    'PHA1QN'    : polarization.xBinnedCountSpectrum,
    'PHA1UN'    : polarization.xBinnedCountSpectrum,
    'CMAP'      : misc.xBinnedMap,
    'MDPMAP'    : polarization.xBinnedMDPMapCube,
    'MDPMAPCUBE': polarization.xBinnedMDPMapCube,
    'PCUBE'     : polarization.xBinnedPolarizationCube,
    'PMAP'      : polarization.xBinnedPolarizationMapCube,
    'PMAPCUBE'  : polarization.xBinnedPolarizationMapCube,
    'PP'        : misc.xBinnedPulseProfile,
    'ARMAP'     : detector.xBinnedAreaRateMap,
    'EFLUX'     : detector.xBinnedAreaEnergyFluxMap,
    'LC'        : misc.xBinnedLightCurve,
    'LTCUBE'    : exposure.xBinnedLivetimeCube
}


def read_binned_file_list(bin_alg, file_list):
    """Read a binned data product file list and return the appropriate data
    structure.
    """
    return BINNING_READ_DICT[bin_alg].from_file_list(file_list)

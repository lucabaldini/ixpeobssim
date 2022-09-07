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

"""This is a variable point source with an un-physical, step-like, time profile.

The range of fluxes explored in different time intervals is appropriate to study
the livetime calculation in xpselect.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.core.spline import xStepFunction
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.time_ import LAUNCH_DATETIME, string_to_met_utc


__model__ = file_path_to_model_name(__file__)

SRC_RA, SRC_DEC = 45., 30.
START_DATE = LAUNCH_DATETIME
START_MET = string_to_met_utc(START_DATE)
DURATION = 50000
STOP_MET = START_MET + DURATION
NUM_INTERVALS = 5
TIME_INTERVALS = numpy.linspace(START_MET, STOP_MET, NUM_INTERVALS + 1)
PL_NORM = numpy.geomspace(0.1, 10., NUM_INTERVALS)
PL_INDEX = 2.
POL_DEG = 0.2
POL_ANG = 40.


pl_norm = xStepFunction(TIME_INTERVALS, PL_NORM, xlabel='MET [s]', ylabel='PL normalization')
spec = power_law(pl_norm, PL_INDEX)
pol_deg = constant(POL_DEG)
pol_ang = constant(numpy.radians(POL_ANG))

src = xPointSource('Point source', SRC_RA, SRC_DEC, spec, pol_deg, pol_ang)
src.set_count_spectrum_params(500, 3, 1)
ROI_MODEL = xROIModel(SRC_RA, SRC_DEC, src)


def display():
    """Display the source model.
    """
    plt.figure('%s PL normalization' % __model__)
    pl_norm.plot(annotate=True)
    setup_gca(grids=True)
    plt.show()



if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

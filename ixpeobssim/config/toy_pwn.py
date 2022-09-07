#!/usr/bin/env python
#
# Copyright (C) 2020, the ixpeobssim team.
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

"""
Toy PWN model---aka known as a gaussian disk with a point source in the middle.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xPointSource, xGaussianDisk, xUniformDisk, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.units_ import arcmin_to_degrees
from ixpeobssim.utils.fmtaxis import fmtaxis


__model__ = file_path_to_model_name(__file__)
ra, dec = 30., 45.

# The pulsar...
pl_norm_pulsar = 1.
pl_index_pulsar = 2.
spec_pulsar = power_law(pl_norm_pulsar, pl_index_pulsar)
pd_pulsar = 0.3
pa_pulsar = 45.
pol_deg_pulsar = constant(pd_pulsar)
pol_ang_pulsar = constant(numpy.radians(pa_pulsar))

#... and the nebula.
pl_norm_nebula = 10.
pl_index_nebula = 2.
radius = arcmin_to_degrees(1.)
spec_nebula = power_law(pl_norm_nebula, pl_index_nebula)
pd_nebula = 0.
pa_nebula = 0.
pol_deg_nebula = constant(pd_nebula)
pol_ang_nebula = constant(numpy.radians(pa_nebula))


pulsar = xPointSource('Pulsar', ra, dec, spec_pulsar, pol_deg_pulsar, pol_ang_pulsar)
nebula = xUniformDisk('Nebula', ra, dec, radius, spec_nebula, pol_deg_nebula, pol_ang_nebula)

ROI_MODEL = xROIModel(ra, dec, pulsar, nebula)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    pass


if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

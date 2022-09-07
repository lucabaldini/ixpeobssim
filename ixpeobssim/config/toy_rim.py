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
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xUniformDisk, xUniformAnnulus, xPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.utils.units_ import arcmin_to_degrees, arcsec_to_degrees
from ixpeobssim.utils.fmtaxis import fmtaxis


__model__ = file_path_to_model_name(__file__)

RA, DEC = 30., 45.


# The big, unpolarized disk.
#r = arcmin_to_degrees(8.)
#spec = power_law(10., 2.)
#pol_deg = constant(0.)
#pol_ang = constant(numpy.radians(0.))
#disk = xUniformDisk('Disk', RA, DEC, r, spec, pol_deg, pol_ang)

# The annulus.
rmin = arcmin_to_degrees(3.5)
rmax = arcmin_to_degrees(3.6)
spec = power_law(5., 2.)
pol_deg = constant(0.5)
pol_ang = constant(numpy.radians(30.))
annulus = xUniformAnnulus('Annulus', RA, DEC, rmin, rmax, spec, pol_deg, pol_ang)

# And a point-like source at the center.
r = arcsec_to_degrees(0.25)
spec = power_law(0.4, 2.)
pol_deg = constant(0.25)
pol_ang = constant(numpy.radians(60.))
source1 = xUniformDisk('Source 1', RA, DEC, r, spec, pol_deg, pol_ang)

# And a point-like source offset.
r = arcsec_to_degrees(0.25)
delta = arcmin_to_degrees(1.95)
spec = power_law(0.1, 2.)
pol_deg = constant(0.25)
pol_ang = constant(numpy.radians(60.))
source2 = xUniformDisk('Source 2', RA + delta, DEC + delta, r, spec, pol_deg, pol_ang)


ROI_MODEL = xROIModel(RA, DEC, annulus, source1, source2)

def display():
    """
    """
    pass

if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

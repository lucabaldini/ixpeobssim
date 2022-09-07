#!/usr/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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

from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.srcmodel.calibsrc import xCalC, xCalibrationROIModel
from ixpeobssim.utils.matplotlib_ import plt, setup_gca


__model__ = file_path_to_model_name(__file__)

rate = 150.

src = xCalC(rate)

ROI_MODEL = xCalibrationROIModel(src)


def display():
    """Display the source model.
    """
    plt.figure('%s morphology' % __model__)
    src.image.plot()
    #plt.plot(energy, spec(energy))
    #setup_gca(xmin=emin, xmax=emax, ymin=spec(emax), logx=True, logy=True,
    #          grids=True, **fmtaxis.spec)


if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

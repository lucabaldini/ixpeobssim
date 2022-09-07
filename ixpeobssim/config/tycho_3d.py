#!/usr/bin/env python
#
# Copyright (C) 2019, the ixpe team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public License as published by
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

import os

from ixpeobssim import IXPEOBSSIM_CONFIG_ASCII, IXPEOBSSIM_CONFIG_FITS
from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.srcmodel.polarization import xStokesSkyCube
from ixpeobssim.srcmodel.roi import xExtendedSource, xROIModel
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.spectrum import load_spectral_spline
from ixpeobssim.srcmodel.img import xFITSImage


__model__ = file_path_to_model_name(__file__)

# Morphology.
ra, dec = 6.340, 64.137
img_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'tycho_4p1_6p1_keV.fits')

# Energy spectrum.
spec_file_path = os.path.join(IXPEOBSSIM_CONFIG_ASCII, 'tycho_total_spectrum.csv')
spec_spline = load_spectral_spline(spec_file_path, delimiter=',', k=1)
spec = lambda E, t: spec_spline(E)

# Polarization
pol_cube = xStokesSkyCube()
inputs = [
    ('polx_0.4_pf_0.30_radial.fits', 'poly_0.4_pf_0.30_radial.fits', 1., None),
    ('polx_0.4_pf_0.30_radial.fits', 'poly_0.4_pf_0.30_radial.fits', 2., 2.83),
    ('polx_0.4_pf_0.60_radial.fits', 'poly_0.4_pf_0.60_radial.fits', 2.83, 4.),
    ('polx_0.4_pf_0.85_radial.fits', 'poly_0.4_pf_0.85_radial.fits', 4.0, 5.66),
    ('polx_0.4_pf_0.90_radial.fits', 'poly_0.4_pf_0.90_radial.fits', 5.66, 8.),
    ('polx_0.4_pf_0.90_radial.fits', 'poly_0.4_pf_0.90_radial.fits', 12., None)
]
for x_file_name, y_file_name, emin, emax in inputs:
    x_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, x_file_name)
    y_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, y_file_name)
    pol_cube.add_layer_xy(x_file_path, y_file_path, emin, emax, rotate=True)
pol_deg = pol_cube.polarization_degree_model()
pol_ang = pol_cube.polarization_angle_model()

# Create the actual ROI model.
tycho = xExtendedSource('Tycho', img_file_path, spec, pol_deg, pol_ang)
ROI_MODEL = xROIModel(ra, dec, tycho)


def display():
    """Display the source model.
    """
    # Energy spectrum
    plt.figure('%s spectrum' % __model__)
    spec_spline.plot()
    setup_gca(xmin=1., xmax=12., logx=True, logy=True, grids=True)

    # Morphology
    plt.figure('%s morphology' % __model__)
    img = xFITSImage(img_file_path)
    img.plot()



if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

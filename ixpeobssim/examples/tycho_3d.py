# Copyright (C) 2016--2022, the ixpeobssim team.
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

import os

import numpy

from ixpeobssim import IXPEOBSSIM_CONFIG_FITS
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
import ixpeobssim.core.pipeline as pipeline
import ixpeobssim.config.tycho_3d as input_model
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color
from ixpeobssim.utils.logging_ import logger


DURATION = 1000000.
MC = False
ANNULI = {'inner': (0., 1.5), 'mid': (1.5, 3.), 'outer': (3., 6.)}


def simulate():
    """Run the simulation and fold the events in phase.
    """
    pipeline.xpobssim(duration=DURATION)


def align():
    """Run the radial alignment tool.

    Note that the center is taken from the WCS in the input polarization models.
    """
    ra0 = 6.331282985
    dec0 = 64.14042889
    file_list = pipeline.file_list()
    pipeline.xpstokesalign(*file_list, ra=ra0, dec=dec0, mode='TAN')
    if not MC:
        return
    # Run the tool using the Monte Carlo position information.
    pipeline.xpstokesalign(*file_list, ra=ra0, dec=dec0, mc=True,
                           suffix='mc', mode='TAN')


def align_from_files():
    """Test the phi alignment from model files.
    """
    file_list = pipeline.file_list()
    x_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'polx_0.4_pf_0.90_radial.fits')
    y_file_path = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'poly_0.4_pf_0.90_radial.fits')
    model_files = [x_file_path, y_file_path]
    pipeline.xpstokesalign(*file_list, mode='XY', modelfiles=model_files)
    if not MC:
        return
    pipeline.xpstokesalign(*file_list, mode='XY', modelfiles=model_files, mc=True,
                           suffix='mc')


def analyze_annuli():
    """Select the external annulus of the Tycho image.
    """
    ra0 = input_model.ra
    dec0 = input_model.dec
    file_list = pipeline.file_list('phialign')
    for label, (rmin, rmax) in ANNULI.items():
        kwargs = dict(ra=ra0, dec=dec0, rad=rmax, innerrad=rmin, suffix=label)
        _list = pipeline.xpselect(*file_list, **kwargs)
        pipeline.xpbin(*_list, algorithm='CMAP')
        pipeline.xpbin(*_list, algorithm='PCUBE', ebins=3, ebinalg='LOG')
    if not MC:
        return
    # Run the analysis using the Monte Carlo position information.
    file_list = pipeline.file_list('phialign', 'mc')
    for label, (rmin, rmax) in ANNULI.items():
        kwargs = dict(ra=ra0, dec=dec0, rad=rmax, innerrad=rmin, suffix=label,
                      mc=True)
        _list = pipeline.xpselect(*file_list, **kwargs)
        pipeline.xpbin(*_list, algorithm='CMAP')
        pipeline.xpbin(*_list, algorithm='PCUBE', ebins=3, ebinalg='LOG')


def analyze_reference_spots():
    """Select and analyze four circular spots on the four sides of tycho to
    control the polarization degree and angle.
    """
    reference_points = {
        'west': (6.16, 64.14),
        'north': (6.34, 64.21),
        'east': (6.51, 64.14),
        'south': (6.34, 64.05)
    }
    file_list = pipeline.file_list()
    for label, (ra, dec) in reference_points.items():
        kwargs = dict(ra=ra, dec=dec, rad=1., suffix=labelprint(pol_deg))
        _list = pipeline.xpselect(*file_list, **kwargs)
        pipeline.xpbin(*_list, algorithm='CMAP')
        pipeline.xpbin(*_list, algorithm='PCUBE', ebins=3, ebinalg='LOG')


def plot():
    """Plot some stuff.
    """
    for label, (rmin, rmax) in ANNULI.items():
        # Plot the modulation cubes.
        _list = pipeline.file_list('phialign_%s' % label, 'pcube')
        pcube = xBinnedPolarizationCube.from_file_list(_list)
        pcube.plot_polarization_degree(label=label)
        # Grab the event file (e.g., for DU #1) after the proper region
        # selection and plot the weighted average of the model.
        file_path = pipeline.file_list('phialign_%s' % label)[0]
        event_file = xEventFile(file_path)
        model_average = event_file.pol_deg_weighted_average(input_model.pol_deg)
        fmt = dict(color=last_line_color(), fmt='-', ls='dashed', errorevery=100)
        model_average.plot(**fmt)
    plt.legend()


def run():
    """Run all.
    """
    simulate()
    align()
    analyze_annuli()
    plot()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline('tycho_3d')

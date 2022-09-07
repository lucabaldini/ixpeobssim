#!/usr/bin/env python
#
# Copyright (C) 2019--2022, the ixpeobssim team.
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

from ixpeobssim import IXPEOBSSIM_IRFGEN
from ixpeobssim.core.fitting import power_law_analytical_fit
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.core.spline import xInterpolatedBivariateSplineLinear
from ixpeobssim.irf.ebounds import ENERGY_STEP, ENERGY_MAX
from ixpeobssim.instrument import COMBINED_DU_ID, du_id_is_valid
from ixpeobssim.utils.logging_ import logger


IRFGEN_MMA_DATA = os.path.join(IXPEOBSSIM_IRFGEN, 'data', 'mma')

MMA_AEFF_FILE_NAME = 'mirror_aeff_calib_20210831.txt'
MMA_AEFF_FILE_PATH = os.path.join(IRFGEN_MMA_DATA, MMA_AEFF_FILE_NAME)

MMA_VIGN_FILE_NAME = 'mirror_vign_v1.txt'
MMA_VIGN_FILE_PATH = os.path.join(IRFGEN_MMA_DATA, MMA_VIGN_FILE_NAME)


def effective_area_spline(du_id, extend_lo=False, extend_hi=False):
    """Return a spline with the on-axis MMA effective area as a function of
    the energy for a given DU.

    This is reading a text file that is periodically updated by MSFC.

    Note that prior to the release of version 6 of the IRF, the text files
    included a single effective-area column with the three MMA combined, and
    this function included an option to perform the division by three.
    As of August 2020, the new files include 5 columns---the energy, the
    effective area for the three DUs, and the total (i.e., sum of the three DUs)
    effective area.

    Also, as of version 6 of the response functions, the effective area is
    defined between 0 and 15 keV (as opposed to 1--12 keV), so that we have
    to do something clever between 0--1 keV and 12--15 keV.

    Arguments
    ---------
    file_name : string
        The name of the text file with the input data
    """
    file_path = MMA_AEFF_FILE_PATH
    logger.info('Reading MMA effective area data from %s...' % file_path)
    try:
        # New file format---five columns, with the energy, the effective area
        # for the three DUs and the sum.
        energy, aeff_du1, aeff_du2, aeff_du3, aeff_total = numpy.loadtxt(file_path, unpack=True)
        # In this case we want to make sure that the sum of the three DUs
        # gives the correct total.
        delta = (aeff_total - (aeff_du1 + aeff_du2 + aeff_du3)) / aeff_total
        assert max(abs(delta)) < 0.02
    except ValueError:
        # Fall back to the old file format, with two columns: energy and
        # total effective area.
        energy, aeff_total = numpy.loadtxt(file_path, unpack=True)
        aeff_du1 = aeff_du2 = aeff_du3 = aeff_total / 3.

    if extend_lo:
        # Since after issue #369 the response functions are defined between 0 and
        # 15 keV, we have to do something below 1 keV and above 12 keV, where the
        # actual MMA effective area is tabulated.
        #
        # Between 0 and 1 keV we simply set the effective area to 0. for energy = 0.,
        # and let the spline interpolate linearly below 1 keV. Since we have
        # essentially zero sensitivity below 1 keV, I don't think the details of what
        # we're doing have any impact, here.
        energy = numpy.append(0., energy)
        aeff_du1 = numpy.append(0., aeff_du1)
        aeff_du2 = numpy.append(0., aeff_du2)
        aeff_du3 = numpy.append(0., aeff_du3)

    def extend_aeff(energy, aeff, energy_ext, fit_emin=10., fit_emax=12.):
        """Small nested function to extend the effective area at high energy.

        Here essentially we fit the high-energy part of the effective area with
        a power low and we extrapolate, imposing continuity in the last
        tabulated point.
        """
        mask = numpy.logical_and(energy >= fit_emin, energy <= fit_emax)
        x = energy[mask]
        y = aeff[mask]
        _, index = power_law_analytical_fit(x, y)
        norm = y[-1] * x[-1]**-index
        return norm * energy_ext**index

    if extend_hi:
        # Above 12 keV we fit the high-energy parte of the MMA effective area with
        # a power law, and use the fit to extrapolate.
        energy_ext = numpy.arange(energy.max() + ENERGY_STEP, ENERGY_MAX, ENERGY_STEP)
        aeff_du1 = numpy.append(aeff_du1, extend_aeff(energy, aeff_du1, energy_ext))
        aeff_du2 = numpy.append(aeff_du2, extend_aeff(energy, aeff_du2, energy_ext))
        aeff_du3 = numpy.append(aeff_du3, extend_aeff(energy, aeff_du3, energy_ext))
        energy = numpy.append(energy, energy_ext)

    # Make sure the du_id passed as an argument is valid.
    assert du_id_is_valid(du_id)
    aeff_dict = {1: aeff_du1, 2: aeff_du2, 3: aeff_du3, COMBINED_DU_ID: aeff_total}
    aeff = aeff_dict[du_id]
    fmt = dict(xlabel='Energy [keV]', ylabel='Effective area [cm$^2$]')
    return xInterpolatedUnivariateSplineLinear(energy, aeff, **fmt)


def effective_area(energy, du_id):
    """Return the on-axis mirror effective area evaluated on a given grid of
    energy points.
    """
    aeff = effective_area_spline(du_id)(energy)
    assert (aeff > 0).all()
    return aeff


def vignetting_spline():
    """Return the vignetting bivariate spline.

    The format of the input file is unfortunate in that the off-axis angles
    are written in the first line, which is commented, and cannote be readily
    parsed with numpy.loadtxt(). We resort instead to opening the file
    "by hand" and parsing the line. The mechanism is fragile, and the assert
    in the code is meant to break things if any change in the file is not
    backward compatible.
    """
    file_path = MMA_VIGN_FILE_PATH
    logger.info('Reading MMA vignetting data from %s...' % file_path)
    data = numpy.loadtxt(file_path)
    energy = data.T[0, :]
    with open(file_path) as f:
        theta = numpy.array([float(val) for val in f.readline().split()[2:]])
    assert numpy.array_equal(theta, numpy.arange(0., 9., 0.5))
    vignetting = data[:,1:]
    fmt = dict(xlabel='Energy [keV]', ylabel='Off-axis angle [arcmin]',
               zlabel='Vignetting')
    return xInterpolatedBivariateSplineLinear(energy, theta, vignetting, **fmt)


def naive_vignetting(energy, theta):
    """Simple, energy-independent, parametrization of the MMA aberration as a
    function of energy and off-axis angle.

    See http://bigfoot.iaps.inaf.it:8080/xwiki/wiki/ixpeglobal/view/Main/IXPEsensitivityFiles/

    Arguments
    ---------
    energy : float or array
       the energy in keV

    theta : float or array (same as energy)
       the off-axis angle in arcmin
    """
    return 1. - 0.0178 * theta - 0.00371 * (theta**2.)


def naive_aberration(energy, theta):
    """Simple, energy-independent, parametrization of the MMA aberration as a
    function of energy and off-axis angle.

    See http://bigfoot.iaps.inaf.it:8080/xwiki/wiki/ixpeglobal/view/Main/IXPEsensitivityFiles/

    Arguments
    ---------
    energy : float or array
       the energy in keV

    theta : float or array (same as energy)
       the off-axis angle in arcmin

    Returns
    -------
    The rms blur radius sigma in arcseconds.
    """
    return 0.0123 * theta + 0.112 * (theta**2.)

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

from __future__ import print_function, division

import os
import numbers

from astropy.io import fits
import numpy

from ixpeobssim.core.spline import xInterpolatedBivariateSpline
from ixpeobssim.core.spline import xUnivariateSpline, xInterpolatedUnivariateSpline
from ixpeobssim import IXPEOBSSIM_IRFGEN
from ixpeobssim.core.modeling import xLine
from ixpeobssim.core.fitting import fit
from ixpeobssim.utils.matplotlib_ import plt, last_line_color, setup_gca
from ixpeobssim.utils.logging_ import logger


# pylint: disable=invalid-name



class xModulationFactorScaling(xInterpolatedBivariateSpline):

    """Convenience class designed to hold the scaling of the modulation fraction
    with energy and pressure, as predicted by the Geant 4 Monte Carlo
    simulation.

    The basic idea is that we run a fine grid of ixpesim simulations for lines
    at different enegries and pressure value, feed them through the standard
    spot analysis and store the measured modulation factor (and associated)
    error into a text file.

    This class is essentially a subclass of a bivariate spline reading the
    file and re-arranging the information in order to make it easy to retrieve
    the modulation factor as a function of pressure and energy.
    """

    def __init__(self, file_name='ixpesim_modf_scan_20200615.txt', correct=True,
                 interactive=False):
        """Constructor.
        """
        file_path = os.path.join(IXPEOBSSIM_IRFGEN, 'data', 'gpd', file_name)
        E, p, modf = self._process_ixpesim_data(file_path, interactive)
        if correct:
            E, p, modf = self._adjust_to_calibration_data(E, p, modf, interactive)
        fmt = dict(xlabel='Energy [keV]', ylabel='Pressure [mbar]',
                   zlabel='Modulation factor', kx=3, ky=1)
        # Impose that the modulation factor for zero energy is zero, no matter
        # what the pressure is. This was added in response to issue #369.
        E = numpy.append(0., E)
        modf = numpy.vstack((numpy.zeros(p.shape), modf))
        xInterpolatedBivariateSpline.__init__(self, E, p, modf, **fmt)

    def value(self, E, p, default=0.45):
        """Alternative to the native __call__() method, to be used at the IRF
        generation stage.

        This is a slightly modified version of the native __call__() method
        where we return an adjustable default value if the energy exceeds the
        maximum energy in the underlying array used for the spline creation.
        Physically, this is used to convey some fictional representation of the
        modulation factor above the Cu K edge, where the analysis needs
        essentially to be written from scratch.

        Note that p needs to be a scalar for this call to work.
        """
        assert isinstance(p, numbers.Number)
        if isinstance(E, numbers.Number):
            E = numpy.array([E])
        mask = E <= self.x.max()
        modf = numpy.zeros(E.shape)
        modf[mask] = xInterpolatedBivariateSpline.__call__(self, E[mask], p)
        modf[numpy.logical_not(mask)] = default
        return modf

    @staticmethod
    def _process_ixpesim_data(file_path, interactive=False):
        """Small routing to process the data into the ixpesim text output.

        Since we tipically run 100k events for each setting in ixpesim, the
        data are just too noisy to make an interpolated bivariate spline
        directly without any pre-processing. Since between 600 and 800 mbar
        the modulation factor as a function of the pressure scales approximately
        linerly at any given energy, we first fit each energy slice with a
        straight, and use the model to set the interpolation points for the
        spline.

        Curiously enough, ixpesim predicts that the slope of the modulation
        factor with energy, in the standard 2--8 keV energy range, is pretty
        much independent from the energy, and is of the order of -1.5e-4 mbar^{-1}.
        That is to say that if the pressure decrease by 100 mbar, the
        modulation factor increases by an additive 1.5% over the IXPE energy band.
        In relative terms, the effect is just less than 5% at 3 keV, and about
        half of that at 6 keV.
        """
        logger.info('Loading modulation factor scan from %s...', file_path)
        energy, pressure, modf, sigma_modf, modf_nocut, sigma_modf_nocut =\
            numpy.loadtxt(file_path, unpack=True)
        logger.info('Done.')
        # Retrieve the underlying energy and pressure vectors.
        energy_grid = numpy.unique(energy)
        pressure_grid = numpy.unique(pressure)
        logger.info('Fitting energy slices...')
        fit_models = []
        # Loop over the energy lines and fit the modulation factor as a function
        # of the pressure with a straight line.
        if interactive:
            plt.figure('Modulation factor fitting')
            setup_gca(xlabel='Pressure [mbar]', ylabel='Modulation factor',
                      xmin=590., xmax=810., ymin=0., ymax=0.7, grids=True)
        for E in energy_grid:
            mask = energy == E
            p = pressure[mask]
            m = modf[mask]
            dm = sigma_modf[mask]
            model = fit(xLine(), p, m, sigma=dm)
            fit_models.append(model)
            if interactive:
                plt.errorbar(p, m, dm, fmt='o', label='%.2f keV' % E)
                plt.plot(p, model(p), color=last_line_color())
                plt.legend()
        # Fit statistics.
        chisq = sum([m.chisq for m in fit_models])
        ndof = sum([m.ndof for m in fit_models])
        logger.info('Overall chisquare: %.2f / %d dof', chisq, ndof)
        # Post-process the list of fit models to get the relevant parameters.
        slope = numpy.array([m.parameter_value('Slope') for m in fit_models])
        sigma_slope = numpy.array([m.parameter_error('Slope') for m in fit_models])
        # Create the ingredients for the bivariate spline for the parent class.
        z = numpy.zeros((len(energy_grid), len(pressure_grid)))
        for i, E in enumerate(energy_grid):
            for j, p in enumerate(pressure_grid):
                z[i, j] = fit_models[i](p)
        return energy_grid, pressure_grid, z

    @staticmethod
    def _adjust_to_calibration_data(E, p, modf, interactive=False):
        """Post-process the ixpesim data going into the final spline to match
        (on average) the calibration data.

        This is essentially loading the calibration data for the flight DUs,
        interpolating the Monte Carlo at the appropriate pressure for the
        measurement time, and calculating an energy-dependent correction
        factor that is in turn applied to the ixpesim data themselves.
        """
        # Select the DU and corresponding calibration pressures.
        du_id = [2, 3]
        calibration_pressure = [725, 680]
        # Build a temporary spline for interpolating the ixpesim values.
        # Note this would be the final spline, if the correction was not applied,
        # and there might be a way of avoiding creating this intermediate
        # product.
        mc_spline = xInterpolatedBivariateSpline(E, p, modf)
        # Load the calibration data and build the correction.
        ratio = numpy.zeros(7)
        for du, pressure in zip(du_id, calibration_pressure):
            file_name = 'modf_du%d_20200605.txt' % du
            file_path = os.path.join(IXPEOBSSIM_IRFGEN, 'data', 'gpd', file_name)
            energy, mu, sigma_mu = numpy.loadtxt(file_path, unpack=True)
            ratio += mu / mc_spline.hslice(pressure)(energy)
        ratio /= len(du_id)
        correction = xInterpolatedUnivariateSpline(energy, ratio, k=1, ext='const')
        if interactive:
            plt.figure('Calibration data correction factor')
            plt.plot(energy, ratio, 'o')
            x = numpy.linspace(1., 9., 250)
            plt.plot(x, correction(x))
            setup_gca(xlabel='Energy [keV]', ylabel='Correction factor', grids=True)
        # Scale the modulation factor values.
        for i, _E in enumerate(E):
            modf[i, :] *= correction(_E)
        return E, p, modf

#!/usr/bin/env python
#
# Copyright (C) 2017--2019, the ixpeobssim team.
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
from ixpeobssim.irfgen.xcom import load_xsection_data
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.core.spline import xInterpolatedPiecewiseUnivariateSpline
from ixpeobssim.core.spline import xInterpolatedBivariateSpline
from ixpeobssim.core.fitting import fit
from ixpeobssim.core.modeling import xLine
from ixpeobssim.irfgen.auxiliary import load_qeff_table, load_allx_rmf_hist,\
    load_lines_rmf_hist, load_pressure_scan_table, AUX_VERSION, AUX_REFERENCE_PRESSURE,\
    lines_qeff_file_path, pscan_modf_file_path
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color
from ixpeobssim.irfgen.constants import BE_DENSITY, AL_DENSITY, dme_density

#pylint: disable=invalid-name


IRFGEN_GPD_DATA = os.path.join(IXPEOBSSIM_IRFGEN, 'data', 'gpd')


def _gpd_data_path(file_name):
    """
    """
    return os.path.join(IRFGEN_GPD_DATA, file_name)


# And now all the relevant nominal detector parameters
WINDOW_BE_THICKNESS = 5.0e-3 # cm (i.e., 50 um)
WINDOW_AL_THICKNESS = 5.3e-6 # cm (i.e., 53 nm)
ABSORPTION_GAP_THICKNESS = 1.0 # cm
GPD_FILL_TEMPERATURE = 20. # degrees C
GPD_FILL_PRESSURE = 800. # mbar
GPD_TYPICAL_ASYMTPTOTIC_PRESSURE = 640. # mbar


# Dictionary of all the contaminants in the Be windows, as per the Materion
# data sheet. The key is the chemical formula of the contaminant, and the value
# if the amount of contaminant, in ppm, relative to the nominal Be density.
# See https://materion.com/-/media/files/electrofusion/eqf3003_purityspecs.pdf
# and https://materion.com/-/media/files/electrofusion/eqf3004_transmittance.pdf
BE_CONTAMINANTS_SPECS = {'Al' : 500.,
                         'BeO': 8000.,
                         'B'  : 3.,
                         'Cd' : 2.,
                         'Ca' : 100.,
                         'C'  : 600.,
                         'Cr' : 100.,
                         'Co' : 10.,
                         'Cu' : 100.,
                         'Fe' : 800.,
                         'Pb' : 20.,
                         'Li' : 3.,
                         'Mg' : 490.,
                         'Mn' : 100.,
                         'Mo' : 20.,
                         'Ni' : 200.,
                         'Si' : 400.,
                         'Ag' : 10.,
                         'Th' : 2.,
                         'U'  : 140.
                         }

# And this is the actual list and concentration of contaminants from the
# certification that came from Materion along with the flight batch of Be foils.
BE_CONTAMINANTS_CERT = {'Al' : 300.,
                        'BeO': 6000.,
                        'B'  : 3.,
                        'Cd' : 2.,
                        'Ca' : 100.,
                        'C'  : 400.,
                        'Cr' : 100.,
                        'Co' : 10.,
                        'Cu' : 100.,
                        'Fe' : 800.,
                        'Pb' : 20.,
                        'Li' : 3.,
                        'Mg' : 490.,
                        'Mn' : 100.,
                        'Mo' : 20.,
                        'Ni' : 100.,
                        'Si' : 300.,
                        'Ag' : 10.,
                        'Th' : 2.,
                        'U'  : 140.
                        }

DEFAULT_BE_CONTAMINANTS = BE_CONTAMINANTS_CERT



def window_be_transparency(energy, contaminants=DEFAULT_BE_CONTAMINANTS,
                           thickness=WINDOW_BE_THICKNESS):
    """Return the transparency of the 50 um Be window.

    Arguments
    ---------
    energy : array-like
        The array of energy values where we calculate the transparency.

    contaminants : dict, optional
        A dictionary containing the list and concentration of contaminants in the
        Be window.
    """
    # Load the Be cross sections and calculate the transparency.
    table = load_xsection_data('Be')
    trans = table.transparency(energy, thickness, BE_DENSITY)
    # If we do not include the contaminant, we're good to go.
    if contaminants is None:
        return trans
    # If we include the contaminants, we need to loop over the proper
    # dictionary and multiply the relative transparency for each element and/or
    # compound.
    for identifier, ppm in contaminants.items():
        table = load_xsection_data(identifier)
        trans *= table.transparency(energy, thickness, BE_DENSITY * ppm * 1.e-6)
    # Ship it!
    return trans


def window_al_transparency(energy, thickness=WINDOW_AL_THICKNESS):
    """Return the transparency of the 53 nm alumination of the Be window.
    """
    table = load_xsection_data('Al')
    return table.transparency(energy, thickness, AL_DENSITY)


def window_transparency(energy, contaminants=DEFAULT_BE_CONTAMINANTS):
    """Return the overall window transparency.
    """
    return window_be_transparency(energy, contaminants) * \
        window_al_transparency(energy)


def photoabsorption_efficiency(energy, temperature, pressure,
                               thickness=ABSORPTION_GAP_THICKNESS):
    """Return the photoabsorption efficiency of the GPD absorption gap.
    """
    table = load_xsection_data('DME')
    density = dme_density(temperature, pressure)
    return table.photoabsorption_efficiency(energy, thickness, density)


def load_ixpesim_ancillary_data(file_name='ixpesim_stdlines_20191109.txt'):
    """Load the ancillary data calculated with the full Geant 4 Monte Carlo
    simulation to inform the generation of the instrument response functions.

    The columns represent, in order:

    * Energy [keV];
    * Trigger efficiency;
    * Probability of extraction of a photoelectron from the Be window;
    * Probability of extraction of a photoelectron from the GEM Cu upper layer;
    * Simulated value of the modulation factor.
    """
    file_path = _gpd_data_path(file_name)
    return numpy.loadtxt(file_path, unpack=True)


def dme_photoemission_frac_spline(pressure=800.):
    """Return a spline representing the fraction of photoelectons extracted
    from the DME as a function of energy.

    This is achieved by calculating the probability of extracting a photoelectron
    from the Be window and the GEM using the Geant 4 Monte Carlo simulation.

    .. warning:
       This function is obsolete and is only kept for diagnostic purposes, use
       the new xQeffDataInterface insterface, instead.
    """
    energy, trg_eff, win_prob, gem_prob, _ = load_ixpesim_ancillary_data()
    dme_prob = 1. - win_prob - gem_prob
    # Implement a first-order correction for the pressure.
    dme_eff = dme_prob * trg_eff * pressure / 800.
    win_eff = win_prob * trg_eff
    gem_eff = gem_prob * trg_eff
    dme_prob = dme_eff / (dme_eff + win_eff + gem_eff)
    return xInterpolatedUnivariateSplineLinear(energy, dme_prob)



class xQeffDataInterface:

    """Basic interface to the post-processed ``ixpesim`` output that is
    necessary to create the arf response functions.

    This is essentially reading the FITS file with the quantum efficiency tabulated
    as a function of the energy for the different conversion types (window, DME,
    and GEM) and caching the necessary arrays in a form that is convenient for
    later use.

    In addition, the class is parametrizing the extraction probabilities for the
    window and the GEM, so that we can readily scale all the relevant quantities
    at different internal pressure in a self-consistent fashion without the
    need to re-run the original simulations.

    .. warning::

       When we updated this to allow for weights, it became clear that the
       efficiency of the active gas volume could not be simply computed from the
       analystic formula (which was the original approach).

       We therefore decided to cache the (energy-dependent) weighting efficiency
       along with the other quantities, based on the consideration that
       its dinamic range is comparatively small, and it is therefore more easily
       amenable to interpolation.
    """

    AL_LINE_ENERGY = 1.557
    CU_LINE_ENERGY = 8.945

    def __init__(self, weight_name=None, aux_version=AUX_VERSION):
        """Constructor.
        """
        self.weight_name = weight_name
        self.aux_version = aux_version
        self.data = load_qeff_table(lines_qeff_file_path(weight_name, aux_version))
        energy = self.data['ENERGY']
        args = energy, self.data['QEFF_WIN'], self.AL_LINE_ENERGY
        self._extr_prob_win_spline = xInterpolatedPiecewiseUnivariateSpline(*args, k=2)
        # In order to model the extraction probability from the GEM we need
        # to deconcolve the effect of the window and the gas absorption, and
        # we need to do that in the same exact conditions of the ixpesim simulation,
        # i.e., no contaminants, 20 degrees and AUX_REFERENCE_PRESSURE mbar.
        # Note that the extraction probability from the GEM, calculated in this
        # fashion, is independent on how we model the Be window, as long as the
        # modeling is consistent between the calculation and the original ixpesim
        # simulations.
        win_trans = window_transparency(energy, contaminants=None)
        dme_qeff = photoabsorption_efficiency(energy, 20., AUX_REFERENCE_PRESSURE)
        extr_prob_gem = self.data['QEFF_GEM'] / (1. - dme_qeff) / win_trans
        args = energy, extr_prob_gem, self.CU_LINE_ENERGY
        self._extr_prob_gem_spline = xInterpolatedPiecewiseUnivariateSpline(*args, k=3)
        # And, in order to support weights, we need to keep track of the
        # weighting efficiency as a function of the energy.
        if weight_name is not None:
            unweighted_data = load_qeff_table(lines_qeff_file_path(None, aux_version))
            ratio_win = self.data['QEFF_WIN'] / unweighted_data['QEFF_WIN']
            ratio_dme = self.data['QEFF_DME'] / unweighted_data['QEFF_DME']
            ratio_gem = self.data['QEFF_GEM'] / unweighted_data['QEFF_GEM']
        else:
            ratio_win = ratio_dme = ratio_gem = numpy.full(energy.shape, 1.)
        self._weight_eff_win_spline = xInterpolatedUnivariateSpline(energy, ratio_win)
        self._weight_eff_dme_spline = xInterpolatedUnivariateSpline(energy, ratio_dme)
        self._weight_eff_gem_spline = xInterpolatedUnivariateSpline(energy, ratio_gem)

    def weight_eff_win(self, energy):
        """Return the window weighting efficiency at a given array of energies.
        """
        return self._weight_eff_win_spline(energy)

    def weight_eff_dme(self, energy):
        """Return the DME weighting efficiency at a given array of energies.
        """
        return self._weight_eff_dme_spline(energy)

    def weight_eff_gem(self, energy):
        """Return the GEM weighting efficiency at a given array of energies.
        """
        return self._weight_eff_gem_spline(energy)

    def window_extraction_prob(self, energy):
        """Return the photoelectron extraction efficiency for the Be window at a
        given energy.

        Mind we clip the output of the spline interpolation to the [0., 1.]
        interval, as the spline might accidentally go below zero at low energy
        due to the noise.
        """
        return numpy.clip(self._extr_prob_win_spline(energy), 0., 1.)

    def window_quantum_efficiency(self, energy):
        """Another name for the same thing :-)

        The window sees all the photons we generate with ixpesim, so the
        denominator is 1, and the quantum efficiency is exactly equal to the
        extraction probability.
        """
        return self.window_extraction_prob(energy)

    def window_absorption_prob(self, energy, temperature=20., pressure=800.,
                               contaminants=DEFAULT_BE_CONTAMINANTS):
        """Return the probabilty for a photon of being absorbed in the window
        at a given energy (or array of energies).
        """
        args = energy, temperature, pressure, contaminants
        return self.window_quantum_efficiency(energy) / self.quantum_efficiency(*args)

    def dme_quantum_efficiency(self, energy, temperature=20., pressure=800.,
                               contaminants=DEFAULT_BE_CONTAMINANTS):
        """Return the quantum efficiency of the DME.
        """
        return window_transparency(energy, contaminants) *\
            photoabsorption_efficiency(energy, temperature, pressure) *\
            self.weight_eff_dme(energy)

    def dme_absorption_prob(self, energy, temperature=20., pressure=800.,
                            contaminants=DEFAULT_BE_CONTAMINANTS):
        """Return the probabilty for a photon of being absorbed in the DME
        at a given energy (or array of energies).
        """
        args = energy, temperature, pressure, contaminants
        return self.dme_quantum_efficiency(*args) / self.quantum_efficiency(*args)

    def gem_extraction_prob(self, energy):
        """Return the photoelectron extraction probability for the GEM at a given energy.

        Mind we clip the output of the spline interpolation to the [0, inf]
        interval, as the spline might accidentally go below zero at low energy
        due to the noise.
        """
        return numpy.clip(self._extr_prob_gem_spline(energy), 0., numpy.inf)

    def gem_quantum_efficiency(self, energy, temperature=20., pressure=800.,
                               contaminants=DEFAULT_BE_CONTAMINANTS):
        """Return the quantum efficiency of the GEM.

        Note that, while the GEM extraction probability is pressure-independent,
        this depends on the pressure due to the effect of the absorptions in the
        gas, changing the number of X-ray photons reaching the GEM top.
        """
        return window_transparency(energy, contaminants) *\
            (1. - photoabsorption_efficiency(energy, temperature, pressure)) *\
            self.gem_extraction_prob(energy)

    def gem_absorption_prob(self, energy, temperature=20., pressure=800.,
                            contaminants=DEFAULT_BE_CONTAMINANTS):
        """Return the probabilty for a photon of being absorbed in the GEM
        at a given energy (or array of energies).
        """
        args = energy, temperature, pressure, contaminants
        return self.gem_quantum_efficiency(*args) / self.quantum_efficiency(*args)

    def quantum_efficiency(self, energy, temperature=20., pressure=800.,
                           contaminants=DEFAULT_BE_CONTAMINANTS):
        """Return the ovearall GPD quantum efficiency, including all the
        absoprtion types.
        """
        args = energy, temperature, pressure
        return self.window_quantum_efficiency(energy) +\
            self.dme_quantum_efficiency(*args, contaminants) +\
            self.gem_quantum_efficiency(*args, contaminants)

    def plot_quantum_efficiency(self, temperature=20., pressure=800.):
        """Plot the quantum efficiency as a function of the energy, disaggregated
        by conversion type.
        """
        energy = numpy.linspace(1., 12., 250)
        args = energy, temperature, pressure
        plt.plot(energy, self.window_quantum_efficiency(energy), label='Window')
        plt.plot(energy, self.dme_quantum_efficiency(*args), label='DME')
        plt.plot(energy, self.gem_quantum_efficiency(*args), label='GEM')
        plt.plot(energy, self.quantum_efficiency(*args), label='All')
        setup_gca(xlabel='Energy [keV]', ylabel='Quantum efficiency @ %d mbar' % pressure,
                  legend=True, grids=True, logy=True, ymin=1.e-4, xmin=1.)

    def plot_extraction_probability(self):
        """Plot the probability of extracting a photoelectron from either the
        window or the GEM as a function of energy.
        """
        energy = numpy.linspace(1., 12., 250)
        plt.plot(energy, self.window_extraction_prob(energy), label='Window')
        plt.plot(energy, self.gem_extraction_prob(energy), label='GEM')
        setup_gca(xlabel='Energy [keV]', ylabel='Photoelectron extraction probability',
                  grids=True, legend=True, xmin=1., ymax=0.01)



def quantum_efficiency(energy, temperature, pressure, weight_name=None, aux_version=AUX_VERSION):
    """Return the overall quantum efficiency of the GPD.

    Arguments
    ---------
    energy : array-like
        The energy array where the quantum efficiency is calculated

    temperature : float
        The GPD filling temperature

    pressure : float
        The GPD gas pressure

    weight_name : str
        The label for the weights to be used.

    aux_version : int
        The version of the auxiliary data products used to model the passive
        conversions.
    """
    qeff_data = xQeffDataInterface(weight_name, aux_version)
    return qeff_data.quantum_efficiency(energy, temperature, pressure)



class xEdispDataInterface:

    """Basic interface to the post-processed ``ixpesim`` output that is
    necessary to create the response matrix.
    """

    KEYS = ('win', 'dme', 'gem')

    def __init__(self, weight_name=None, aux_version=AUX_VERSION):
        """Constructor.
        """
        # Warning---should the weight name be set for the qeff interface?
        self.qeff_data = xQeffDataInterface(weight_name, aux_version)
        self.win_hist = load_allx_rmf_hist('win', weight_name, aux_version)
        self.dme_hist = load_lines_rmf_hist('dme', weight_name, aux_version)
        self.gem_hist = load_allx_rmf_hist('gem', weight_name, aux_version)

    def combine(self, temperature=20., pressure=800., contaminants=DEFAULT_BE_CONTAMINANTS):
        """Combine the window, DME and GEM response functions in the proper
        proportions (depending on the target temperature and pressure, as well
        as the window contaminants) to create the actual response function.
        """
        hist = self.dme_hist.empty_copy()
        energy = hist.bin_centers(1)
        args = energy, temperature, pressure, contaminants
        win_prob = self.qeff_data.window_absorption_prob(*args)
        dme_prob = self.qeff_data.dme_absorption_prob(*args)
        gem_prob = self.qeff_data.gem_absorption_prob(*args)
        data = numpy.zeros(hist.shape)
        for i, E in enumerate(energy):
            row = self.win_hist.content[:, i] * win_prob[i] +\
              self.dme_hist.content[:, i] * dme_prob[i] +\
              self.gem_hist.content[:, i] * gem_prob[i]
            row /= row.sum()
            data[:, i] = row
        hist.set_content(data)
        return hist



class xModfDataInterface(xInterpolatedBivariateSpline):

    """Basic interface to the post-processed ``ixpesim`` output parametrizing
    the modulatuion factor as a function of the pressure.
    """

    def __init__(self, weight_name=None, aux_version=AUX_VERSION):
        """
        """
        file_path = pscan_modf_file_path(weight_name, aux_version)
        self.pressure, self.energy, self.mu, self.mu_err = load_pressure_scan_table(file_path)
        fitted_mu = self._fit_data()
        fmt = dict(xlabel='Energy [keV]', ylabel='Pressure [mbar]',
                   zlabel='Modulation factor', kx=2, ky=1)
        xInterpolatedBivariateSpline.__init__(self, self.energy, self.pressure, fitted_mu, **fmt)

    def _fit_data(self, interactive=False):
        """Fit the underlying pressure scan data.

        Here we plot the modulation factor as a function of the pressure for
        each value of the energy and fit the points with a straight line, then
        build a numpy array of the same shape of `self.mu` where the fit
        model is evaluated at each pressure.
        """
        fitted_mu = numpy.zeros(self.mu.shape)
        for i, E in enumerate(self.energy):
            y = self.mu[i]
            dy = self.mu_err[i]
            model = fit(xLine(), self.pressure, y, sigma=dy)
            fitted_mu[i] = model(self.pressure)
            if interactive:
                plt.figure('Modulation factor @ %.2f keV' % E)
                plt.errorbar(self.pressure, y, dy, fmt='o')
                model.plot(color=last_line_color())
                model.stat_box()
                setup_gca(xlabel='Pressure [mbar]', ylabel='Modulation factor', grids=True)
        return fitted_mu

    def pressure_slice(self, pressure):
        """Return a slice at the given pressure, in the form of a piecewise
        interpolated spline.
        """
        args = self.energy, self.hslice(pressure)(self.energy), xQeffDataInterface.CU_LINE_ENERGY
        return xInterpolatedPiecewiseUnivariateSpline(*args, k=2)



def modulation_factor(energy, pressure, weight_name=None, aux_version=AUX_VERSION):
    """Return the modulation factor calculated on a grid of points.
    """
    return xModfDataInterface(weight_name, aux_version).pressure_slice(pressure)(energy)

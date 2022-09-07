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

"""Library for producing auxiliary files for the IRF generation.
"""

from __future__ import print_function, division

import os
import sys

from astropy.io import fits
import numpy
import scipy
from scipy.optimize import curve_fit

from ixpeobssim import IXPEOBSSIM_IRFGEN
from ixpeobssim.core.fitting import fit_histogram, fit_gaussian_iterative
from ixpeobssim.core.hist import xHistogram1d, xHistogram2d
from ixpeobssim.core.modeling import xFitModelBase, xLogNormal, xGeneralizedGaussian, xHat
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.irf.ebounds import channel_to_energy, PI_BINNING
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca

if sys.flags.interactive:
    plt.ion()


#pylint: disable=invalid-name, too-many-locals, too-many-arguments


AUX_OUTPUT_FOLDER = os.path.join(IXPEOBSSIM_IRFGEN, 'data', 'gpd')
AUX_VERSION = 3
AUX_REFERENCE_PRESSURE = 700.
AUX_REFERENCE_TEMPERATURE = 20.
AUX_WINDOW_CONTAMINANTS = None
AUX_CUSTOM_SPEC_FILE_NAME = 'allx_spec_v001.dat'
AUX_ABSORPTION_LABELS = ('win', 'dme', 'gem')
AUX_WEIGHT_NAMES = (None, 'alpha075')


class _xEdispModelLogNormal(xFitModelBase):

    """Energy dispersion fitting model.

    This is the combination of a log-normal distribution for the main peak, and
    a hat-shaped bridge for the tail, with the possible addition of a low-energy
    bump.
    """

    #pylint: disable=arguments-differ

    PARAMETER_NAMES = ('Normalization', 'PeakFrac', 'Peak', 'Sigma', 'Shape',
                       'Endpoint', 'BumpFrac')
    PARAMETER_DEFAULT_VALUES = (1., 1., 100., 20., 0.05, 20., 0.01)
    PARAMETER_DEFAULT_BOUNDS = ([0., 0.6, 10., 3., 0.001, 2., 0.],
                                [numpy.inf, 1., 350., 50., 1., 50., 0.1])
    DEFAULT_PLOTTING_RANGE = (0., 375.)
    BUMP_SIGMA = 5.

    @staticmethod
    def value(x, norm, peak_frac, peak, sigma, shape, endpoint, bump_frac):
        """Overloaded value() method.
        """
        return norm * (peak_frac * xLogNormal.value(x, 1., peak, sigma, shape) +\
            (1 - peak_frac - bump_frac) * xHat.value(x, 1., endpoint, 3., peak, sigma) +\
            bump_frac * scipy.stats.norm.pdf(x, endpoint, xEdispModelLogNormal.BUMP_SIGMA))


class xEdispModelLogNormal(xFitModelBase):

    """Energy dispersion fitting model.

    This is the combination of a log-normal distribution for the main peak, and
    a hat-shaped bridge for the tail, with the possible addition of a low-energy
    bump.
    """

    #pylint: disable=arguments-differ

    PARAMETER_NAMES = ('Normalization', 'PeakFrac', 'Peak', 'Sigma', 'Shape',
                       'Endpoint', 'BumpFrac', 'BumpSigma')
    PARAMETER_DEFAULT_VALUES = (1., 1., 100., 20., 0.05, 20., 0.01, 3.)
    PARAMETER_DEFAULT_BOUNDS = ([0., 0.6, 10., 3., 0.001, 2., 0., 3.],
                                [numpy.inf, 1., 350., 50., 1., 50., 0.1, 8.])
    DEFAULT_PLOTTING_RANGE = (0., 375.)

    @staticmethod
    def value(x, norm, peak_frac, peak, sigma, shape, endpoint, bump_frac, bump_sigma):
        """Overloaded value() method.
        """
        return norm * (peak_frac * xLogNormal.value(x, 1., peak, sigma, shape) +\
            (1 - peak_frac - bump_frac) * xHat.value(x, 1., endpoint, 3., peak, sigma) +\
            bump_frac * scipy.stats.norm.pdf(x, endpoint, bump_sigma))



class xEdispModelGenGauss(xFitModelBase):

    """Energy dispersion fitting model.

    This is the combination of a generalized Gaussian distribution for the main peak, and
    a hat-shaped bridge for the tail, with the possible addition of a low-energy
    bump.
    """

    #pylint: disable=arguments-differ

    PARAMETER_NAMES = ('Normalization', 'PeakFrac', 'Peak', 'Sigma', 'Beta', 'BumpFrac')
    PARAMETER_DEFAULT_VALUES = (1., 0.8, 100., 20., 1., 0.2)
    PARAMETER_DEFAULT_BOUNDS = ([0., 0.25, 15., 2., 1., 0.],
                                [numpy.inf, 1., 350., 50., 3., 0.6])
    DEFAULT_PLOTTING_RANGE = (0., 375.)
    BUMP = 18.
    BUMP_SIGMA = 6.

    @staticmethod
    def value(x, norm, peak_frac, peak, sigma, beta, bump_frac):
        """Overloaded value() method.
        """
        bump = xEdispModelGenGauss.BUMP
        bump_sigma = xEdispModelGenGauss.BUMP_SIGMA
        return norm * (peak_frac * xGeneralizedGaussian.value(x, 1., peak, sigma, beta) +\
            bump_frac * scipy.stats.norm.pdf(x, bump, bump_sigma) +\
            (1. - peak_frac - bump_frac) * xHat.value(x, 1., bump, bump_sigma, peak, sigma))



def _aux_file_path_base(base_name, label=None, version=AUX_VERSION):
    """Return the path to the output file to write the histogram data for a given
    absorption type and aux file version.

    This is the main book-keeping function for file paths related to data products
    for auxiliary files.
    """
    if label is None:
        file_name = 'aux_%s_v%03d.fits' % (base_name, version)
    else:
        file_name = 'aux_%s_%s_v%03d.fits' % (base_name, label, version)
    return os.path.join(AUX_OUTPUT_FOLDER, file_name)


def _edisp_file_path(base_name, abs_label, weight_name, aux_version=AUX_VERSION):
    """Base function for edisp pre-processed files.
    """
    label = abs_label
    if weight_name is not None:
        label = '%s_%s' % (label, weight_name)
    return _aux_file_path_base(base_name, label, aux_version)


def allx_edisp_file_path(abs_label, weight_name, aux_version=AUX_VERSION):
    """Return the path to the allx edisp data file path for a given absorption
    type and set of weights.
    """
    return _edisp_file_path('allx_edisp', abs_label, weight_name, aux_version)


def allx_qeff_file_path(weight_name, aux_version=AUX_VERSION):
    """Return the path to the allx lines data file path for a given absorption
    type and set of weights.
    """
    return _aux_file_path_base('allx_qeff', weight_name, aux_version)


def allx_pha_model_file_path(aux_version=AUX_VERSION):
    """Return the path to the file with the PHA model.
    """
    return _aux_file_path_base('allx_pha_model', None, aux_version)


def allx_rmf_file_path(abs_label, weight_name=None, aux_version=AUX_VERSION):
    """Return the path to the FITS file(s) to be used to compose the response matrix.
    """
    return _edisp_file_path('allx_rmf', abs_label, weight_name, aux_version)


def lines_edisp_file_path(abs_label, weight_name, aux_version=AUX_VERSION):
    """Return the path to the lines edisp data file path for a given absorption
    type and set of weights.
    """
    return _edisp_file_path('lines_edisp', abs_label, weight_name, aux_version)


def lines_qeff_file_path(weight_name, aux_version=AUX_VERSION):
    """Return the path to the lines qeff data file path for a given absorption
    type and set of weights.
    """
    return _aux_file_path_base('lines_qeff', weight_name, aux_version)


def lines_rmf_file_path(abs_label, weight_name=None, aux_version=AUX_VERSION):
    """Return the path to the FITS file(s) to be used to compose the response matrix.
    """
    return _edisp_file_path('lines_rmf', abs_label, weight_name, aux_version)


def pscan_modf_file_path(weight_name, aux_version=AUX_VERSION):
    """Return the path to the modulation factor data from the analysis of the
    pressure scan.
    """
    return _aux_file_path_base('pscan_modf', weight_name, aux_version)


def calculate_efficiency(n, N):
    """Calculate the efficiency, given a number of trials and success, and
    attach the classical binomial error to it.
    """
    eff = n / N
    eff_err = numpy.sqrt(n * (N - n) / N**3.)
    return eff, eff_err


def _process_headers(hdu_list):
    """Run over all the headers for a FITS file and collect all the relevant
    keywords (i.e., starting with APID for the primary header, or those coming
    after the EXTNAME key for the other extensions) into a single dictionary.
    """
    info = {}
    for i, hdu in enumerate(hdu_list):
        keys = tuple(hdu.header.keys())
        if i == 0:
            index = keys.index('APID')
        else:
            index = keys.index('EXTNAME') + 1
        keys = keys[index:]
        info.update({key: hdu.header.get(key) for key in keys})
    # Note we pop out the HISTORY keyword.
    info.pop('HISTORY')
    return info


def absorption_type_masks(mc_absz, min_dme_absz=0.832, max_dme_absz=10.830):
    """Return a set of three masks for the window, DME and GEM absorption.

    This is based on the Monte Carlo z coordinate of the absorption point, and
    it is not the most elegant solution, as the correctness of the mask relies
    on the fact that we use the right bounds for the thickness of the absorption
    gap used in the simulations. (Since the latter is essentially never changed,
    this is less horrible than it might seems at a first glance.)

    Moving forward, fixing
    https://bitbucket.org/ixpesw/gpdsw/issues/222
    and use the exact information might be the way to go.
    """
    window_mask = mc_absz > max_dme_absz
    gem_mask = mc_absz < min_dme_absz
    dme_mask = numpy.logical_and(mc_absz >= min_dme_absz, mc_absz <= max_dme_absz)
    return window_mask, dme_mask, gem_mask


def _fit_pi_peak(hist, true_energy):
    """Small fitting routine to find the peak of an histogram for the purpose of
    converting energy to pulse invariant in ixpesim photon lists.
    """
    # First Gaussian fit to gauge the mean and standard deviation of the
    # histogram.
    p0 = (hist.num_entries(), 3300., 300.)
    gaussian = fit_gaussian_iterative(hist, p0=p0, num_sigma_left=1., num_sigma_right=2.)
    # Refined, log-normal fit. Mind the cut on the left of the peak is
    # quite aggressive, but was necessary to make the fit stable.
    # Also, above the Cu edge, where there is a lot of garbage on the left of the
    # main peak, we need to initialize the shape of the log-normal to a small
    # value for the fit to converge.
    if true_energy < 8.945:
        p0 = (hist.num_entries(), gaussian.Peak, gaussian.Sigma, 0.05)
    else:
        p0 = (hist.num_entries(), gaussian.Peak, gaussian.Sigma, 1.e-5)
    xmin = gaussian.Peak - 1.5 * gaussian.Sigma
    try:
        return fit_histogram(xLogNormal(), hist, p0=p0, xmin=xmin)
    except RuntimeError:
        logger.info('Minimum value for the log-normal fit: %.1f', xmin)
        logger.info('Starting parameter values for the log-normal fit: %s', p0)
        plt.figure('Fitting debug @ %.2f keV' % true_energy)
        hist.plot()
        gaussian.plot()
        plt.show()
        logger.warning('Falling back to Gaussian approximation')
        model = xLogNormal()
        model.set_parameters(*p0)
        return model


def _norm_pha_hist1d(mc_energy, pha):
    """Create a one-dimensional histogram of the normalized PHA.
    """
    binning = numpy.linspace(0., 5000., 201)
    norm_pha = pha / mc_energy
    return xHistogram1d(binning, xlabel='PHA / Energy [keV$^{-1}$]').fill(norm_pha)


def _norm_pha_hist2d(mc_energy, pha, energy_binning=None):
    """Create a two-dimensional histogram of the normalized PHA.
    """
    pha_binning = numpy.linspace(0., 5000., 201)
    if energy_binning is None:
        energy_binning = numpy.linspace(mc_energy.min(), mc_energy.max(), 50)
    fmt = dict(xlabel='PHA / True energy [keV$^{-1}$]', ylabel='True energy [keV]')
    norm_pha = pha / mc_energy
    return xHistogram2d(pha_binning, energy_binning, **fmt).fill(norm_pha, mc_energy)


def calculate_pi_broadband(mc_energy, pha, energy_binning=None):
    """Calculate the event pulse invariant for a broadband simulation.

    This is trying to account for the bulk of the non linearities by
    fitting the normalized PHA in different energy slices and then fitting
    the scale factor as a function of the energy with an erf function.
    (The effect is of the order of 2% from 2 keV to 12 keV, at about
    3000 ADC counts per keV nominal.)
    """
    # Create a 2-dimensional histogram of the true energy vs. PHA.
    hist2d = _norm_pha_hist2d(mc_energy, pha, energy_binning)
    energy = hist2d.vslice(0).bin_centers()
    # Fit all the horizontal slices with the baseline model.
    fit_models = [_fit_pi_peak(h, E) for h, E in zip(hist2d.hslices(), energy)]
    scale = numpy.array([model.Mean for model in fit_models])
    scale_err = numpy.array([model.Sigma for model in fit_models])

    def _pha_fit_function(x, const, norm, mean, sigma):
        """Fit function for the conversion scale between energy and PHA.

        Note this is purely phenomenological.
        """
        #pylint: disable=no-member
        return const + norm * scipy.special.erf((x - mean) / sigma)

    # Note we only use the energies below the copper edge to fit the energy
    # dependence of the PHA scaling factor.
    mask = energy < 8.9
    _x = energy[mask]
    _y = scale[mask]
    _dy = scale_err[mask]
    p0 = (3000., 100., 5., 10.)
    popt, _ = curve_fit(_pha_fit_function, _x, _y, sigma=_dy, p0=p0)
    # Debugging plots.
    plt.figure('PI conversion matrix')
    hist2d.plot()
    plt.figure('PI conversion model')
    plt.plot(energy, scale, 'o')
    plt.plot(energy, _pha_fit_function(energy, *popt))
    setup_gca(xlabel='Energy [keV]', ylabel='PI scale', grids=True)
    # Create the model to convert PHA to PI.
    _E = channel_to_energy(PI_BINNING)
    _pha = _E * _pha_fit_function(_E, *popt)
    pha_model = xInterpolatedUnivariateSpline(_pha, PI_BINNING, xlabel='PHA', ylabel='PI')
    pi = pha_model(pha)
    rec_energy = channel_to_energy(pi)
    return rec_energy, pi, pha_model


def calculate_pi_line(pha, aux_version=AUX_VERSION):
    """Convert the pha column into pulse invariant for a line simulation.

    Note that, for line simulations, we use the PHA model derived from the
    allx data set.
    """
    pha_model = load_pha_model(aux_version)
    pi = pha_model(pha)
    rec_energy = channel_to_energy(pi)
    return rec_energy, pi, pha_model


def calculate_pi(mc_energy, pha):
    """Perform the conversion from PHA to PI.

    This is necessary because the PHA values in the ixpesim simulation are
    expressed in the terms of the native detector ADC counts, and the
    scale factor to the actual energy channels in the response functions is
    arbitrary. In real life the conversion will be based on the ground calibration
    with monochromatic sources, and this is meant to be a sensible proxy
    for that process for operating on Monte Carlo data sets.

    Note this operates differently depending on whether we're working with
    line data or a continoum spectrum.

    Args
    ----
    mc_energy : array_like
        The array of (true) energies in keV

    pha : array_like
        The array of measured PHA values.

    Returns
    -------
    rec_energy : array_like
        The array of reconstructed energies (in keV) corresponding to the input PHAs.

    pi : array_like
        The array of pulse invariants corresponding to the input PHAs.

    pha_model : xInterpolatedUnivariateSpline object
        The model to operate the conversion between PHAs and PIs.
    """
    # If the dynamics in energy is less than 10%, then we'll call it a line.
    if mc_energy.max() / mc_energy.min() < 1.1:
        return calculate_pi_line(pha)
    return calculate_pi_broadband(mc_energy, pha)


def _calculate_alpha(event_data):
    """Calculate alpha, i.e., the event ellipticity, given a series of events.

    Args
    ----
    event_data : dict
        Dictionary containing the columns of an input photon list.
    """
    m2l, m2t = event_data['TRK_M2L'], event_data['TRK_M2T']
    return (m2l - m2t) / (m2l + m2t)


def event_weights(event_data, label='alpha075'):
    """Book-keeping function mapping text labels to weights.
    """
    if label is None:
        return None
    if label == 'alpha075':
        return _calculate_alpha(event_data)**0.75
    logger.error('Weight label "%s" not cought', label)
    return None


def load_event_data(*file_list, absorption_masks=True, pi=True):
    """Load the data from the (reconstructed) event files generates with ixpesim.

    This is chaining together all the data in the input files---for allx simulations
    you should pass all the files to a single function call, while for lines
    you should be open them one by one.

    The event data are returned in the form of a dictionary indexed by the
    corresponding column names in the EVENTS and MONTE_CARLO extensions of the
    photon lists. Metadata are stored into a separate info dictionary.
    """
    # pylint: disable=no-member
    logger.info('Loading event data...')
    mc_cols = ('ENERGY', 'ABS_Z', 'PE_PHI')
    rec_cols = ('PHA', 'PHI1', 'PHI2', 'TRK_M2L', 'TRK_M2T')
    _num_gen = 0
    _num_evt = 0
    data = {}
    # Loop over the files and retrieve the actual vectors
    for i, file_path in enumerate(file_list):
        logger.info('Opening input file %s...', file_path)
        with fits.open(file_path) as hdu_list:
            mc = hdu_list['MONTE_CARLO'].data
            evt = hdu_list['EVENTS'].data
            _num_gen += hdu_list['MONTE_CARLO'].header['NUM_GEN']
            _num_evt += hdu_list['MONTE_CARLO'].header['NUM_EVT']
            if i == 0:
                # For the first file we have to store the metadata and populate
                # the output dictionary...
                info = _process_headers(hdu_list)
                for col_name in mc_cols:
                    data[col_name] = mc[col_name]
                for col_name in rec_cols:
                    data[col_name] = evt[col_name]
            else:
                # ...while for all the remaining files we keep appending to the
                # original arrays.
                for col_name in mc_cols:
                    data[col_name] = numpy.append(data[col_name], mc[col_name])
                for col_name in rec_cols:
                    data[col_name] = numpy.append(data[col_name], evt[col_name])
    logger.info('Done, %d event(s) read out / %d generated.', _num_evt, _num_gen)
    # Update the info dict.
    info['NUM_GEN'] = _num_gen
    info['NUM_EVT'] = _num_evt
    if absorption_masks:
        # Calculate the absorption masks and add them to the data dict.
        data['WIN_MASK'], data['DME_MASK'], data['GEM_MASK'] =\
            absorption_type_masks(data['ABS_Z'])
    if pi:
        # Perform the PI conversion and update the data dict.
        data['REC_ENERGY'], data['PI'], data['PHA_MODEL'] =\
            calculate_pi(data['ENERGY'], data['PHA'])
    # Calculate the necessary weights.
    weight_dict = {label: event_weights(data, label) for label in AUX_WEIGHT_NAMES}
    return data, info, weight_dict


def write_pha_model(file_path, model, info):
    """Write the PHA to PI conversion model table to file.
    """
    # pylint: disable=no-member
    logger.info('Writing PHA to PI conversion model to %s...', file_path)
    col_names = ('PHA', 'PI')
    data = {'PHA': model.x, 'PI': model.y}
    cols = [fits.Column(name=col, array=data[col], format='E') for col in col_names]
    table = fits.BinTableHDU.from_columns(cols)
    table.name = 'PHA_MODEL'
    hdu_list = fits.HDUList([fits.PrimaryHDU(), table])
    for key, value in info.items():
        hdu_list[0].header.set(key, value)
    hdu_list.writeto(file_path, overwrite=True)
    logger.info('Done.')


def load_pha_model(aux_version=AUX_VERSION):
    """Load a PHA to PI model from file.
    """
    file_path = allx_pha_model_file_path(aux_version)
    # pylint: disable=no-member
    logger.info('Loading PHA to PI conversion model from %s...', file_path)
    with fits.open(file_path) as hdu_list:
        data = hdu_list['PHA_MODEL'].data
        pha = data['PHA']
        pi = data['PI']
    return xInterpolatedUnivariateSpline(pha, pi, xlabel='PHA', ylabel='PI')


def load_allx_edisp_data(abs_label, weight_name=None, aux_version=AUX_VERSION):
    """Load the (raw) response matrix data created with mkauxallx.py.
    """
    file_path = allx_edisp_file_path(abs_label, weight_name, aux_version)
    return xHistogram2d.from_file(file_path)


def load_lines_edisp_data(abs_label, weight_name=None, aux_version=AUX_VERSION):
    """Load the quantum efficiency data from file.
    """
    file_path = lines_edisp_file_path(abs_label, weight_name, aux_version)
    logger.info('Loading edisp data from %s...', file_path)
    with fits.open(file_path) as hdu_list:
        data = hdu_list['EDISP_HISTS'].data
    energy = data['ENERGY']
    hists = []
    for i, _ in enumerate(energy):
        hist = xHistogram1d(PI_BINNING, xlabel='Pulse invariant')
        entries = data['DATA'][i]
        hist.set_content(entries, entries)
        hists.append(hist)
    return energy, hists


def write_qeff_table(file_path, data, info):
    """Write a quantum efficiency table to file.
    """
    # pylint: disable=no-member
    logger.info('Writing quantum efficiency data to %s...', file_path)
    col_names = ('ENERGY', 'QEFF', 'QEFF_ERR', 'QEFF_WIN', 'QEFF_WIN_ERR',
                 'QEFF_DME', 'QEFF_DME_ERR', 'QEFF_GEM', 'QEFF_GEM_ERR')
    cols = [fits.Column(name=col, array=data[col], format='E') for col in col_names]
    table = fits.BinTableHDU.from_columns(cols)
    table.name = 'QEFF'
    hdu_list = fits.HDUList([fits.PrimaryHDU(), table])
    for key, value in info.items():
        hdu_list[0].header.set(key, value)
    hdu_list.writeto(file_path, overwrite=True)
    logger.info('Done.')


def load_qeff_table(file_path):
    """Load the quantum efficiency data from file.

    This is reading the binary table stored in the FITS files and returning
    the columns in the form of a Python dictionary indexed by column name.
    """
    # pylint: disable=no-member
    logger.info('Loading quantum efficiency binned data from %s...', file_path)
    with fits.open(file_path) as hdu_list:
        data = hdu_list['QEFF']
        data = {col.name: data.data[col.name] for col in data.columns}
    return data


def write_pressure_scan_table(file_path, data, info):
    """Write a pressure table to file.
    """
    # pylint: disable=no-member
    logger.info('Writing pressure scan data to %s...', file_path)
    col_names = ('ENERGY', 'PRESSURE', 'MU1', 'MU1_ERR', 'MU2', 'MU2_ERR')
    cols = [fits.Column(name=col, array=data[col], format='E') for col in col_names]
    table = fits.BinTableHDU.from_columns(cols)
    table.name = 'PSCAN'
    hdu_list = fits.HDUList([fits.PrimaryHDU(), table])
    for key, value in info.items():
        hdu_list[0].header.set(key, value)
    hdu_list.writeto(file_path, overwrite=True)
    logger.info('Done.')


def load_pressure_scan_table(file_path, col_name='MU2'):
    """Load the pressure scan data from file.

    This is reading the binary table stored in the FITS files and arranging them
    in an order that is suitable for interpolation in pressure and energy.
    """
    # pylint: disable=no-member
    logger.info('Loading pressure scan binned data from %s...', file_path)
    with fits.open(file_path) as hdu_list:
        data = hdu_list['PSCAN']
        data = {col.name: data.data[col.name] for col in data.columns}

    def _filter_colum(p):
        """Small nested function to retrieve the modulation factor as a
        function of the energy at a given pressure.

        (Note the output values are sorted by ascending energy.)
        """
        mask = data['PRESSURE'] == p
        E = data['ENERGY'][mask]
        idx = numpy.argsort(E)
        val = data[col_name][mask][idx]
        err = data['%s_ERR' % col_name][mask][idx]
        return val, err

    pressure_grid = numpy.unique(data['PRESSURE'])
    energy_grid = numpy.unique(data['ENERGY'])
    mu = numpy.array([], dtype=float)
    mu_err = numpy.array([], dtype=float)
    for p in pressure_grid:
        val, err = _filter_colum(p)
        mu = numpy.append(mu, val)
        mu_err = numpy.append(mu_err, err)
    shape_ = (len(pressure_grid), len(energy_grid))
    mu = mu.reshape(shape_).T
    mu_err = mu_err.reshape(shape_).T
    return pressure_grid, energy_grid, mu, mu_err


def load_allx_rmf_hist(abs_label, weight_name=None, aux_version=AUX_VERSION):
    """Load the post-processed histogram with the rmf data.
    """
    return xHistogram2d.from_file(allx_rmf_file_path(abs_label, weight_name, aux_version))


def load_lines_rmf_hist(abs_label, weight_name=None, aux_version=AUX_VERSION):
    """Load the post-processed histogram with the rmf data.
    """
    return xHistogram2d.from_file(lines_rmf_file_path(abs_label, weight_name, aux_version))

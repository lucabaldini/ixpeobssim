#!/usr/bin/env python
#
# Copyright (C) 2018--2020, the ixpeobssim team.
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

from astropy.io import fits

import xspec

from ixpeobssim import IXPEOBSSIM_XSPEC
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.irf.caldb import irf_folder_path
from ixpeobssim.utils.os_ import check_input_file
from ixpeobssim.binning.base import peek_binning_algorithm, xBinnedFileBase
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, residual_plot, xStatBox
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.utils.matplotlib_ import last_line_color, nlog_errorbars


# pylint: disable=invalid-name

STAT_METHODS = ['chi', 'cstat', 'lstat', 'pgstat', 'pstat', 'whittle']
POLARIMETRIC_ADDITIVE_MODELS = ['apolconst', 'apollin', 'apolpow']


class xXspecSpectrumManager(dict):

    """Small facility to keep track of the spectra being loaded into XSPEC.

    This is essentially a dict where the type of the binned spectra and the
    DU Ids are mapped into the data group in the XSPEC memory, so that we
    can easily retrieve the plot data at any time.
    """

    def register(self, file_path, data_group):
        """Register a spectrum.

        Note we peek into the binned file to get an hold on the binning algorithm
        and the DU ID.
        """
        input_file = xBinnedFileBase(file_path)
        bin_alg = input_file.binning_algorithm()
        du_id = input_file.du_id()
        args = os.path.basename(file_path), bin_alg, du_id, data_group
        logger.info('Registering [...]%s (%s, DU %d) at data group %d', *args)
        self[self._key(bin_alg, du_id)] = data_group

    @staticmethod
    def _key(bin_alg, du_id):
        """Build the indexing key from the binning algorithm and the DU ID.
        """
        return (bin_alg, du_id)

    def data_group(self, bin_alg, du_id, default=None):
        """Return the XSPEC data group mapping to a particular binning algorithm
        and DU ID.
        """
        return self.get(self._key(bin_alg, du_id), default)

    def spectrum_types(self):
        """Return a tuple with all the binning algorithms containing at
        least a spectrum.
        """
        return tuple(set([key[0] for key in self]))



# Global variable for the book-keeping of the data sets loaded in memory.
# This is cleared up and re-populated at each load_input_files() call.
_SPECTRUM_MANAGER = xXspecSpectrumManager()



class xXspecPlotData:

    """Small container class storing the XSPEC plot data.

    Plot data are typically retrieved by function calls to the global xspec.Plot
    object, and this contained provides a convenient interface to cache plot
    data for later use. The container encapsulates the x and y values, the
    errors on the y values and, if available, the model calculated in
    correspondence of the x array.
    """

    def __init__(self, energy, data, errors, model=None):
        """Constructor.
        """
        self.energy = numpy.array(energy)
        self.data = numpy.array(data)
        self.errors = numpy.array(errors)
        if model is not None:
            self.model = numpy.array(model)
        else:
            self.model = None

    def __len__(self):
        """Return the lenght of the underlying data vector.
        """
        return len(self.energy)

    def save(self, file_path):
        """To be implemented.
        """
        raise NotImplementedError

    def __str__(self):
        """String representation.
        """
        if self.model is not None:
            text = ''
        else:
            text = ' not'
        return '%d point(s) found, model%s available' % (len(self), text)



class xXspecFitData:

    """Small containter class to store the output of an XSPEC spectral fit.

    This includes the parameter names, best-fit values and errors, the test
    statistics and the number of degrees of freedom.
    """

    def __init__(self, model):
        """Constructor.
        """
        self.test_statistic = xspec.Fit.testStatistic
        self.stat_test = xspec.Fit.statTest
        self.dof = xspec.Fit.dof
        self.__index = 0
        self.model_expression = model.expression
        self.par_names = []
        self.__par_index_dict = {}
        self.par_values = []
        self.par_errors = []
        self.par_low_bounds = []
        self.par_high_bounds = []
        self.par_status_codes = []
        for i in range(model.nParameters):
            par = model(i + 1)
            self.par_names.append(par.name)
            self.__par_index_dict[par.name] = i
            self.par_values.append(par.values[0])
            self.par_errors.append(par.sigma)
            low_bound, high_bound, status_code = par.error
            self.par_low_bounds.append(low_bound)
            self.par_high_bounds.append(high_bound)
            self.par_status_codes.append(status_code)

    def __call__(self, identifier):
        """Overloaded method returning the parameter value and error at a given
        index.
        """
        if isinstance(identifier, str):
            assert identifier in self.par_names
            identifier = self.__par_index_dict[identifier]
        return self.par_values[identifier], self.par_errors[identifier]

    def par_value(self, identifier):
        """Return the parameter value at a given index.
        """
        if isinstance(identifier, str):
            assert identifier in self.par_names
            identifier = self.__par_index_dict[identifier]
        return self.par_values[identifier]

    def par_error(self, identifier):
        """Return the parameter error at a given index.
        """
        if isinstance(identifier, str):
            assert identifier in self.par_names
            identifier = self.__par_index_dict[identifier]
        return self.par_errors[identifier]

    def __len__(self):
        """Return the number of parameters in the model.
        """
        return len(self.par_names)

    def __iter__(self):
        """Implementation of the iteration protocol.
        """
        return self

    def __next__(self):
        """Implementation of the iteration protocol.
        """
        if self.__index < len(self):
            i = self.__index
            self.__index += 1
            return self.par_names[i], self.par_values[i], self.par_errors[i],\
                   self.par_low_bounds[i], self.par_high_bounds[i], self.par_status_codes[i]
        # I don't know that this is a clever thing to do, here---need to
        # study the Python iteration protocol better.
        self.__index = 0
        raise StopIteration()

    def next(self):
        """Alternative method for Python 2 compatibility.
        """
        return self.__next__()

    def __str__(self):
        """String formatting.
        """
        text = 'Fit model: %s (%s = %.2f / %d)\n' %\
            (self.model_expression, self.stat_test, self.test_statistic, self.dof)
        for name, value, error, low_bound, high_bound, status_code in self:
            if error > 0:
                text += '%15s: %.3e +/- %.3e' % (name, value, error)
                if low_bound != high_bound:
                    pos = high_bound - value
                    neg = value - low_bound
                    text += ' (+%.3e / -%.3e) %s' % (pos, neg, status_code)
                text += '\n'
            else:
                text += '%15s: %.3e     (frozen)\n' % (name, value)
        text.strip('\n')
        return text

    def stat_box(self, **kwargs):
        """Return a stat box that can be overlaid onto a given plot.

        Note that when the error command has been run after the fit, we do
        include the average of the (asymmetric) error bars in the stat box.
        """
        stat_box = xStatBox(**kwargs)
        stat_box.add_entry('Fit model: %s' % self.model_expression)
        stat_box.add_entry('TS (%s) = %.2f / %d dof' %\
                           (self.stat_test, self.test_statistic, self.dof))
        for name, value, error, low_bound, high_bound, status_code in self:
            if low_bound != high_bound:
                error = 0.5 * (high_bound - low_bound)
            stat_box.add_entry(name, value, error)
        return stat_box



def compare_fit_data(fit_data, target_values, threshold=5.):
    """Compare the best-fit parameter values with the input model.

    Arguments
    ---------
    fit_data : xXspecFitData instance
        The fit output.

    target_values : array_like
        The corresponding values from the input model.

    threshold : float
        The threshold (in units of sigma) for determining the agreement
        between the best fit parameters and the model inputs.

    Return
    ------
    Zero if the bes-fit parameters agree with the input, and an error code > 0
    otherwise.
    """
    error_code = 0
    logger.info('Comparing fit output with input model...')
    for i, (name, value, error, _, _, _) in enumerate(fit_data):
        target = target_values[i]
        # This is catching possible frozen parameters.
        if error > 0:
            delta = abs((value - target) /  error)
        else:
            delta = 0
        msg = '%s = %.4e +/- %.4e, target = %.4e (%.2f sigma)' %\
              (name, value, error, target, delta)
        if delta > threshold:
            error_code += 1
            logger.error(msg)
        else:
            logger.info(msg)
    return error_code


def reset():
    """Global reset.
    """
    xspec.AllData.clear()
    # This is needed to avoid XSPEC prompting the user for, e.g., missing
    # response files, which we take care of programmaticaly.
    xspec.Xset.allowPrompting = False


def add_model_string(key, value):
    """Add a key,value pair of strings to XSPEC's internal database.
    (simple wrapper around the xspec.Xset.addModelString() method.)

    This database provides a way to pass string values to certain model functions
    (e.g., POW_EMIN and POW_EMAX) which are hardcoded to search for "key".
    (See the XSPEC manual description for the "xset" command for a table showing
    model/key usage.)

    If the key,value pair already exists, it will be replaced with the new entries.
    """
    logger.info('Setting "%s" to %s in the XSPEC internal database...', key, value)
    xspec.Xset.addModelString(key, str(value))


def add_model_strings(**kwargs):
    """Facility to add multiple model string to xspec.
    """
    for key, value in kwargs.items():
        add_model_string(key, value)


def load_input_files(*file_list):
    """Read a set of PHA1 input files and create the corresponding
    xspec.Spectrum objects.
    """
    # Reset XSPEC and clear the book-keeping dictionary.
    reset()
    _SPECTRUM_MANAGER.clear()
    logger.info('Loading input files...')
    for i, file_path in enumerate(file_list):
        _SPECTRUM_MANAGER.register(file_path, i + 1)
        logger.info('Loading binned count spectrum from %s...', file_path)
        spectrum = xspec.Spectrum(file_path)

        # At this point XSPEC has already looked into the RESPFILE e ANCRFILE
        # keywords in the file header and assigned the corresponing response
        # files to the spectrum. This might not be what the user wants if the
        # paths are absolute paths on a different machine, in which case the
        # spectrum will have no reponse associated and spectrum.response will
        # raise an exception.

        try:
            spectrum.response
        except:
            logger.info('Spectrum response invalid, using local CALDB')

            def _irf_path(file_path, irf_type):
                """Small nested function to load a given response file
                (identified by file name) from the local caldb, independently
                from the absolute path specified in the binned PHA1 file header.
                """
                file_name = os.path.basename(file_path)
                file_path = os.path.join(irf_folder_path(irf_type), file_name)
                check_input_file(file_path)
                logger.info('Using %s file %s...', irf_type, file_path)
                return file_path

            with fits.open(file_path) as hdu_list:
                respfile = hdu_list['SPECTRUM'].header['RESPFILE']
                ancrfile = hdu_list['SPECTRUM'].header['ANCRFILE']
                spectrum.response = _irf_path(respfile, 'rmf')
                if ancrfile.endswith('arf'):
                    spectrum.response.arf = _irf_path(ancrfile, 'arf')
                if ancrfile.endswith('mrf'):
                    spectrum.response.arf = _irf_path(ancrfile, 'mrf')
                if ancrfile.endswith('modf'):
                    spectrum.response.arf = _irf_path(ancrfile, 'modf')

    logger.info('Done, %d file(s) loaded.', len(file_list))


def select_energy_range(emin=None, emax=None):
    """Apply a global mask on the energy value for the channels to be fitted.
    """
    if emin is not None:
        logger.info('Ignoring channels below %.3f keV...', emin)
        xspec.AllData.ignore('**-%f' % emin)
    if emax is not None:
        logger.info('Ignoring channels above %.3f keV...', emax)
        xspec.AllData.ignore('%f-**' % emax)


def load_local_models(package_name='ixpeobssim'):
    """Load all the ixpeobssim local models.

    We introduced this function, and called it automatically, version 12.0.0,
    based on the idea that whoever had XSPEC installed also had the XSPEC
    local models installed. This turned out to be a glaring usability issue,
    as there are many possible ways one could be unable to compile the local
    models, see  https://bitbucket.org/ixpesw/ixpeobssim/issues/350
    As of version 12.1.0 we wrap the XSPEC call in a try / except, and happily
    proceed even if the local models cannot be loaded.
    """
    args = package_name, IXPEOBSSIM_XSPEC
    logger.info('Loading XSPEC local models from the "%s" package in %s...', *args)
    try:
        xspec.AllModels.lmod(*args)
    except:
        logger.error('Could not XSPEC load local models.')
        logger.error('(This might indicate that you need to compile them.)')
        logger.error('See the documentation for more info about XSPEC support.')
    logger.info('Done.')


def setup_fit_model(expression, par_values=None):
    """Create a model for spectral fitting in XSPEC.

    Note we abort if the model cannot be created---this can happen, e.g.,
    if one is trying to use local models without having compiled them. This
    seems as the only sane option, at this point, as the result of moving on
    would depend on the internal XSPEC state (is there any active model?) and
    this would be likely more confusing to the user than stating straight away
    what the problem is.
    """
    expression = expression.strip()
    logger.info('Setting up model "%s"...', expression)
    try:
        model = xspec.Model(expression)
    except:
        # Opsss, cannot create the model...
        logger.error('Could not create model "%s"', expression)
        logger.error('Did you forget to compile the local models?')
        abort('Cannot proceed with fitting')
    if par_values is not None:
        for i, par_value in enumerate(par_values):
            par_index = i + 1
            logger.info('Setting parameter %d to %.2f..', par_index, par_value)
            model(par_index).values = par_value
    # If needed, freeze the phony normalization for additive, purely polarimetric
    # models.
    if expression in POLARIMETRIC_ADDITIVE_MODELS:
        logger.info('Freezing the normalization to 1. for additive model %s...', expression)
        # Note that the normalization is the last parameters, and we start
        # counting from 1.
        parameter = model(model.nParameters)
        parameter.values = 1.
        parameter.frozen = True
    return model


def set_parameter(par_index, par_value, model_id=1):
    """Set the model parameter value corresponding to provided index.
    """
    model = xspec.AllModels(model_id)
    parameter = model(par_index)
    parameter.values = par_value
    return parameter


def fix_parameter(par_index, par_value, model_id=1):
    """Fix a model parameter to a given value.
    """
    parameter = set_parameter(par_index, par_value, model_id)
    parameter.frozen = True


def set_parameter_range(par_index, min_value, max_value, model_id=1):
    """Set the range for a parameter.
    """
    model = xspec.AllModels(model_id)
    parameter = model(par_index)
    values = parameter.values
    values[2] = min_value
    values[3] = min_value
    values[4] = max_value
    values[5] = max_value
    parameter.values = values
    return parameter


def fit(stat_method='chi', max_iterations=10, error=True):
    """Perform a spectral fit.
    """
    assert stat_method in STAT_METHODS
    xspec.Fit.nIterations = max_iterations
    xspec.Fit.statMethod = stat_method
    xspec.Fit.perform()
    if error:
        calculate_confidence_intervals()
    xspec.Fit.show()
    return current_fit_output()


def calculate_confidence_intervals(model_id=1, delta_fit_statistics=1.):
    """Calculate the confidence intervals for a fit.

    This is running the error XSPEC command, see
    https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node78.html
    """
    model = xspec.AllModels(model_id)
    arg_string = '%.5f 1-%d' % (delta_fit_statistics, model.nParameters)
    xspec.Fit.error(arg_string)


def current_fit_output(model_id=1):
    """Return the current fit output.
    """
    model = xspec.AllModels(model_id)
    return xXspecFitData(model)


def retrieve_plot_data(bin_alg, du_id):
    """Retrieve the plot data (and, if available, the corresponding model)
    for a given spectrum.

    Warning
    -------
    Note that the plot data within XSPEC are intented to be normalized wrt the
    energy---i.e., the input from the FITS table is divided by the width of the
    energy bin, and the units get an additional keV^{-1}.
    """
    group = _SPECTRUM_MANAGER.data_group(bin_alg, du_id)
    if group is None:
        return
    xspec.Plot.xAxis = 'keV'
    xspec.Plot('data')
    data = [xspec.Plot.x(group), xspec.Plot.y(group), xspec.Plot.yErr(group)]
    try:
        data.append(xspec.Plot.model(group))
    except:
        pass
    data = xXspecPlotData(*data)
    logger.info('Done, %s', data)
    return data


def plot_normalized_counts(plot_data, du_id, **kwargs):
    """Plot the normalized counts.

    Note that we use full circles for positive values and empty circles for
    negative ones, so that we can plot Stokes parameters (that can go negative)
    in logarithmic scale.
    """
    kwargs.setdefault('zorder', 1)
    kwargs.setdefault('ms', 5)
    kwargs.setdefault('fmt', 'o')
    kwargs['label'] = 'DU %d' % du_id
    nlog_errorbars(plot_data.energy, plot_data.data, plot_data.errors, **kwargs)
    if plot_data.model is not None:
        plt.plot(plot_data.energy, plot_data.model, zorder=2, color='black')


def plot_residuals(plot_data):
    """Plot the model residuals.
    """
    res = (plot_data.data - plot_data.model) / plot_data.errors
    plt.errorbar(plot_data.energy, res, 1., fmt='o', zorder=1, ms=5)
    plt.plot(plot_data.energy, numpy.zeros(plot_data.energy.shape),
             zorder=2, color='black')
    setup_gca(grids=True, **fmtaxis.ene_pulls_sigma)


def _ylabel(bin_alg):
    """Return the units for a given binned spectrum.
    """
    label = 'Normalized %s counts' % bin_alg
    if bin_alg in ('PHA1', 'PHA1Q', 'PHA1U'):
        return '%s [s$^{-1}$ keV$^{-1}$]' % label
    if bin_alg in ('PHA1QN', 'PHA1UN'):
        return '%s [keV$^{-1}$]' % label


def plot(figure_name='XSPEC spectrum', logy=True):
    """Custom, matplotlib-based implementation of the standard XSPEC
    spectral-fit plot.
    """
    for bin_alg in _SPECTRUM_MANAGER.spectrum_types():
        ax1, ax2 = residual_plot('%s %s' % (figure_name, bin_alg))
        for du_id in DU_IDS:
            data = retrieve_plot_data(bin_alg, du_id)
            if data is None:
                continue
            plt.sca(ax1)
            plot_normalized_counts(data, du_id)
            plt.sca(ax2)
            plot_residuals(data)
        plt.sca(ax1)
        setup_gca(ylabel=_ylabel(bin_alg), grids=True, logy=logy)
        fit_data = current_fit_output()
        fit_data.stat_box(position='lower left').plot()
        plt.legend()


def energy_flux(emin=2., emax=8.):
    """Return the value of the energy flux and of the photon flux between
    emin and emax.
    """
    xspec.AllModels.calcFlux('%.3f %.3f err' % (emin, emax))
    return xspec.AllData(1).flux


def sample_spectral_model(expression, parameters, emin=1., emax=12.,
                          num_points=250, name=None, source_id=1):
    """Sample the values for a generic XSPEC model, given an expression and a set
    of parameters.

    Most notably, this is used to feed into ixpeobssim time-independent
    XSPEC spectral models.

    Arguments
    ---------
    expression : str
        The model expression string, using full component names.

    name : str
        The model name.

    source_id : int
        The source number.

    parameters : dict or list
        The model parameters.

    emin, emax : double
        Energy limits, in keV.

    num_points : int
        The number of points to sample the spectrum.
    """
    if name is None:
        name = 'Model%d' % (len(xspec.AllModels.sources) + 1)
    model = xspec.Model(expression, name, source_id, parameters)
    parameter_names = [model(i + 1).name for i in range(model.nParameters)]
    # Set the energy grid for sampling the spectrum.
    # In XSPEC language, this is effectively the same thing as defining a linear
    # binning in energy and then taking the bin centers.
    # Since on our end we are interested in spectra at precise energy points, we
    # add the proper padding, here, so that in the end the output energy array
    # is what we mean with the input parameters---within numerical accuracy.
    padding = 0.5 * (emax - emin) / num_points
    emin -= padding
    emax += padding
    binning = numpy.linspace(emin, emax, num_points + 1)
    xspec.AllModels.setEnergies("%.12f %.12f %i lin" % (emin, emax, num_points))
    xspec.Plot.device = "/null"
    xspec.Plot('model %s' % name)
    energy = numpy.array(xspec.Plot.x())
    # Mind this is the actual flux integrated over the bin and divided by the
    # bin width---not the flux calculated at the bin center, so an
    # a-posteriori correction might be needed.
    flux = numpy.array(xspec.Plot.model())
    xspec.AllModels.setEnergies('reset')
    return binning, energy, flux, parameter_names

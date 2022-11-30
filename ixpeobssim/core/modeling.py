#!/usr/bin/env python
#
# Copyright (C) 2015, the ixpeobssim team.
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


import numpy
import scipy.stats

from ixpeobssim.utils.math_ import fold_angle_deg, fold_angle_rad
from ixpeobssim.utils.matplotlib_ import plt, xStatBox
from ixpeobssim.utils.logging_ import abort


class xFitModelBase:

    """Base class for a fittable model.

    This base class isn't really doing anything useful, the idea being that
    actual models that can be instantiated subclass xFitModelBase overloading
    the relevant class members.

    The class features a number of static members that derived class
    should redefine as needed:

    * ``PARAMETER_NAMES`` is a list containing the names of the model
      parameters. (It goes without saying thet its length should match the
      number of parameters in the model.)
    * ``PARAMETER_DEFAULT_VALUES`` is a list containing the default values
      of the model parameters. When a concrete model object is instantiated
      these are the values being attached to it at creation time.
    * ``PARAMETER_DEFAULT_BOUNDS`` is a tuple containing the default values
      of the parameter bounds to be used for the fitting. The values in the
      tuple are attached to each model object at creation time and are
      intended to be passed as the ``bounds`` argument of the
      ``scipy.optimize.curve_fit()`` function. From the ``scipy`` documentation:
      Lower and upper bounds on independent variables. Defaults to no bounds.
      Each element of the tuple must be either an array with the length equal
      to the number of parameters, or a scalar (in which case the bound is
      taken to be the same for all parameters.) Use np.inf with an appropriate
      sign to disable bounds on all or some parameters. By default
      models have no built-in bounds.
    * ``DEFAULT_PLOTTING_RANGE`` is a two-element list with the default support
      (x-axis range) for the model. This is automatically updated at runtime
      depending on the input data when the model is used in a fit.
    * ``DEFAULT_STAT_BOX_POSITION`` is the default location of the stat box for
      the model, see the :meth:`gpdswpy.matplotlib_.stat_box()` function for
      all the details.

   In addition, each derived class should override the following things:

    * the ``value(x, *args)`` static method: this should return the value of
      the model at a given point for a given set of values of the underlying
      parameters;
    * (optionally) the ``jacobian(x, *args)`` static method. (If defined, this
      is passed to the underlying fit engine allowing to reduce the number of
      function calls in the fit; otherwise the jacobian is calculated
      numerically.)

    Finally, if there is a sensible way to initialize the model parameters
    based on a set of input data, derived classes should overload the
    ``init_parameters(xdata, ydata, sigma)`` method of the base class, as the
    latter is called by fitting routines if no explicit array of initial values
    are passed as an argument. The default behavior of the class method defined
    in the base class is to do nothing.

    See :class:`gpdswpy.modeling.xGaussian` for a working example.
    """

    PARAMETER_NAMES = ()
    PARAMETER_DEFAULT_VALUES = ()
    PARAMETER_DEFAULT_BOUNDS = (-numpy.inf, numpy.inf)
    DEFAULT_PLOTTING_RANGE = (0., 1.)
    DEFAULT_STAT_BOX_POSITION = 'upper right'

    def __init__(self):
        """Constructor.

        Here we initialize the class members holding the best-fit parameter
        values and the associated covariance matrix, see the
        :pymeth:`gpdswpy.modeling.xFitModelBase.reset()` method.

        We also create a (private) look-up table holding the correspondence
        between the parameter names and the corresponding position in the
        parameter list and we cache the default support for the model for
        plotting purposes.
        """
        assert len(self.PARAMETER_NAMES) == len(self.PARAMETER_DEFAULT_VALUES)
        self.__parameter_dict = {}
        for i, name in enumerate(self.PARAMETER_NAMES):
            self.__parameter_dict[name] = i
        self.reset()

    def name(self):
        """Return the model name.
        """
        return self.__class__.__name__

    def reset(self):
        """Reset all the necessary stuff.

        This method initializes all the things that are necessry to keep
        track of a parametric fit.

        * the parameter values are set to what it specified in
          ``PARAMETER_DEFAULT_VALUES``;
        * the covatiance matrix is initialized to a matrix of the proper
          dimension filled with zeroes;
        * the minimum and maximum values of the independent variable
          (relevant for plotting) are set to the values specified in
          ``DEFAULT_PLOTTING_RANGE``;
        * the model bounds are set to the values specified in
          ``PARAMETER_DEFAULT_BOUNDS``;
        * the chisquare and number of degrees of freedom are initialized to
          -1.
        """
        self.parameters = numpy.array(self.PARAMETER_DEFAULT_VALUES, dtype='d')
        self.covariance_matrix = numpy.zeros((len(self), len(self)), dtype='d')
        self.xmin, self.xmax = self.DEFAULT_PLOTTING_RANGE
        self.bounds = self.PARAMETER_DEFAULT_BOUNDS
        self.chisq = -1.
        self.ndof = -1

    def reduced_chisquare(self):
        """Return the reduced chisquare.
        """
        if self.ndof > 0:
            return self.chisq / self.ndof
        return -1.

    def __getattr__(self, name):
        """Short-hand method to retrieve the parameter values by name.

        Note that we manipulate the attribute name by capitalizing the
        first letter and replacing underscores with spaces.
        """
        name = name.title().replace('_', ' ')
        if name in self.PARAMETER_NAMES:
            return self.parameter_value(name)
        else:
            raise AttributeError

    @staticmethod
    def value(x, *parameters):
        """Eval the model at a given point and a given set of parameter values.

        Warning
        -------
        This needs to be overloaded in any derived class for the thing to do
        something sensible.
        """
        raise 'value() not implemented'

    def integral(self, edges):
        """Calculate the integral of the model within pre-defined edges.

        Note that this assumes that the derived class provides a suitable
        ``cdf()`` method.
        """
        try:
            return self.cdf(edges[1:], *self.parameters) - \
                self.cdf(edges[:-1], *self.parameters)
        except Exception as e:
            abort('%s.integral() not implemened (%s)' % (self.name, e))

    def rvs(self, size=1):
        """Return random variates from the model.

        Note that this assumes that the derived class provides a suitable
        ``ppf()`` method.
        """
        try:
            return self.ppf(numpy.random.random(size), *self.parameters)
        except Exception as e:
            abort('%s.rvs() not implemened (%s)' % (self.name(), e))

    def __call__(self, x, *parameters):
        """Return the value of the model at a given point and a given set of
        parameter values.

        Note that unless the proper number of parameters is passed to the
        function call, the model is evaluated at the best-fit parameter values.

        The function is defined with this signature because it is called
        with a set of parameter values during the fit process, while
        tipically we want to evaluate it with the current set of parameter
        values after the fact.
        """
        if len(parameters) == len(self):
            return self.value(x, *parameters)
        else:
            return self.value(x, *self.parameters)

    def __parameter_index(self, name):
        """Convenience method returning the index within the parameter vector
        for a given parameter name.
        """
        assert(name in self.PARAMETER_NAMES)
        return self.__parameter_dict[name]

    def parameter_value(self, name):
        """Return the parameter value by name.
        """
        index = self.__parameter_index(name)
        return self.parameters[index]

    def parameter_error(self, name):
        """Return the parameter error by name.
        """
        index = self.__parameter_index(name)
        return numpy.sqrt(self.covariance_matrix[index][index])

    def parameter_values(self):
        """Return the vector of parameter values.
        """
        return self.parameters

    def parameter_errors(self):
        """Return the vector of parameter errors.
        """
        return numpy.sqrt(self.covariance_matrix.diagonal())

    def parameter_status(self):
        """Return the complete status of the model in the form of a tuple
        of tuples (parameter_name, parameter_value, parameter_error).

        Note this can be overloaded by derived classes if more information
        needs to be added.
        """
        return tuple(zip(self.PARAMETER_NAMES, self.parameter_values(),
                         self.parameter_errors()))

    def set_parameters(self, *pars):
        """Set all the parameter values.

        Note that the arguments must be passed in the right order.
        """
        self.parameters = numpy.array(pars, dtype='d')

    def set_parameter(self, name, value):
        """Set a parameter value.
        """
        index = self.__parameter_index(name)
        self.parameters[index] = value

    def init_parameters(self, xdata, ydata, sigma):
        """Assign a sensible set of values to the model parameters, based
        on a data set to be fitted.

        Note that in the base class the method is not doing anything, but it
        can be reimplemented in derived classes to help make sure the
        fit converges without too much manual intervention.
        """
        pass

    def set_plotting_range(self, xmin, xmax):
        """Set the plotting range.
        """
        self.xmin = xmin
        self.xmax = xmax

    def plot(self, *parameters, **kwargs):
        """Plot the model.

        Note that when this is called with a full set of parameters, the
        self.parameters class member is overwritten so that the right values
        can then be picked up if the stat box is plotted.
        """
        if len(parameters) == len(self):
            self.parameters = parameters
        display_stat_box = kwargs.pop('display_stat_box', False)
        x = numpy.linspace(self.xmin, self.xmax, 1000)
        y = self(x, *parameters)
        plt.plot(x, y, **kwargs)
        if display_stat_box:
            self.stat_box(**kwargs)

    def stat_box(self, position=None, plot=True, **kwargs):
        """Plot a ROOT-style stat box for the model.
        """
        if position is None:
            position = self.DEFAULT_STAT_BOX_POSITION
        box = xStatBox(position)
        box.add_entry('Fit model: %s' % self.name())
        box.add_entry('Chisquare', '%.1f / %d' % (self.chisq, self.ndof))
        for name, value, error in self.parameter_status():
            box.add_entry(name, value, error)
        if plot:
            box.plot(**kwargs)
        return box

    def __len__(self):
        """Return the number of model parameters.
        """
        return len(self.PARAMETER_NAMES)

    def __add__(self, other):
        """Add two models.

        Warning
        -------
        This is highly experimental and guaranteed to contain bugs. Enjoy.
        """
        m1 = self
        m2 = other
        xmin = min(m1.DEFAULT_PLOTTING_RANGE[0], m2.DEFAULT_PLOTTING_RANGE[0])
        xmax = max(m1.DEFAULT_PLOTTING_RANGE[1], m2.DEFAULT_PLOTTING_RANGE[1])
        name = '%s + %s' % (m1.__class__.__name__, m2.__class__.__name__)
        # In order to correctly propagate the parameter boundaries of the two
        # input models to their sum, we need to handle the case where one or
        # both of them do not define their PARAMETER_DEFAULT_BOUNDS, relying
        # on the base class default. In that case we simply repeat those
        # default values for as many parameters as there are in that model.
        lower_bounds = []
        upper_bounds = []
        for m in (m1, m2):
            if m.bounds == xFitModelBase.PARAMETER_DEFAULT_BOUNDS:
                lower_bounds += (len(m.parameters) * [m.bounds[0]])
                upper_bounds += (len(m.parameters) * [m.bounds[1]])
            else:
                lower_bounds += m.bounds[0]
                upper_bounds += m.bounds[1]

        class _model(xFitModelBase):

            PARAMETER_NAMES = m1.PARAMETER_NAMES + m2.PARAMETER_NAMES
            PARAMETER_DEFAULT_VALUES = m1.PARAMETER_DEFAULT_VALUES + \
                                       m2.PARAMETER_DEFAULT_VALUES
            DEFAULT_PLOTTING_RANGE = (xmin, xmax)
            PARAMETER_DEFAULT_BOUNDS = (lower_bounds, upper_bounds)

            def __init__(self):
                self.__class__.__name__ = name
                xFitModelBase.__init__(self)

            @staticmethod
            def value(x, *parameters):
                return m1.value(x, *parameters[:len(m1)]) +\
                    m2.value(x, *parameters[len(m1):])

            @staticmethod
            def jacobian(x, *parameters):
                return numpy.hstack((m1.jacobian(x, *parameters[:len(m1)]),
                                     m2.jacobian(x, *parameters[len(m1):])))

        return _model()

    def __str__(self):
        """String formatting.
        """
        text = '%s model (chisq/ndof = %.2f / %d)' % (self.__class__.__name__,
                                                      self.chisq, self.ndof)
        for name, value, error in self.parameter_status():
            text += '\n%15s: %.5e +- %.5e' % (name, value, error)
        return text


class xConstant(xFitModelBase):

    """Constant model.

    .. math::
      f(x; C) = C
    """

    PARAMETER_NAMES = ('Constant',)
    PARAMETER_DEFAULT_VALUES = (1.,)
    DEFAULT_PLOTTING_RANGE = (0., 1.)

    @staticmethod
    def value(x, constant):
        """Overloaded value() method.
        """
        return numpy.full(x.shape, constant)

    @staticmethod
    def jacobian(x, constant):
        """Overloaded jacobian() method.
        """
        d_constant = numpy.full((len(x),), 1.)
        return numpy.array([d_constant]).transpose()

    def cdf(self, x):
        """Overloaded cdf() method.
        """
        return self.Constant * x

    def ppf(self, q):
        """Overloaded ppf() method.
        """
        return self.xmin + q * (self.xmax - self.xmin)

    def init_parameters(self, xdata, ydata, sigma):
        """Overloaded init_parameters() method.
        """
        self.set_parameter('Constant', numpy.mean(ydata))



class xLine(xFitModelBase):

    """Straight-line model.

    .. math::
      f(x; m, q) = mx + q
    """

    PARAMETER_NAMES = ('Intercept', 'Slope')
    PARAMETER_DEFAULT_VALUES = (1., 1.)
    DEFAULT_PLOTTING_RANGE = (0., 1.)

    @staticmethod
    def value(x, intercept, slope):
        """Overloaded value() method.
        """
        return intercept + slope * x

    @staticmethod
    def jacobian(x, intercept, slope):
        """Overloaded jacobian() method.
        """
        d_intercept = numpy.full((len(x),), 1.)
        d_slope = x
        return numpy.array([d_intercept, d_slope]).transpose()



class xGaussian(xFitModelBase):

    """One-dimensional Gaussian model.

    .. math::
      f(x; A, \\mu, \\sigma) = A e^{-\\frac{(x - \\mu)^2}{2\\sigma^2}}
    """

    PARAMETER_NAMES = ('Amplitude', 'Peak', 'Sigma')
    PARAMETER_DEFAULT_VALUES = (1., 0., 1.)
    PARAMETER_DEFAULT_BOUNDS = ([0., -numpy.inf, 0], [numpy.inf] * 3)
    DEFAULT_PLOTTING_RANGE = (-5., 5.)
    SIGMA_TO_FWHM = 2.3548200450309493

    @staticmethod
    def value(x, amplitude, peak, sigma):
        """Overloaded value() method.
        """
        return amplitude * numpy.exp(-0.5 * ((x - peak)**2. / sigma**2.))

    @staticmethod
    def der_amplitude(x, amplitude, peak, sigma):
        """Return the amplitude derivative of the function, to be used in the
        calculation of the Jacobian.
        """
        return numpy.exp(-0.5 / sigma**2. * (x - peak)**2.)

    @staticmethod
    def der_peak(x, amplitude, d_amplitude, peak, sigma):
        """Return the peak derivative of the function, to be used in the
        calculation of the Jacobian.

        Note that we pass the pre-calculated values of the amplitude derivatives
        in order not to repeat the same calculation more times than strictly
        necessary.
        """
        return amplitude * d_amplitude * (x - peak) / sigma**2.

    @staticmethod
    def der_sigma(x, amplitude, d_amplitude, peak, sigma):
        """Return the sigma derivative of the function, to be used in the
        calculation of the Jacobian.

        Note that we pass the pre-calculated values of the amplitude derivatives
        in order not to repeat the same calculation more times than strictly
        necessary.
        """
        return amplitude * d_amplitude * (x - peak)**2. / sigma**3.

    @staticmethod
    def jacobian(x, amplitude, peak, sigma):
        """Overloaded jacobian() method.
        """
        d_amplitude = xGaussian.der_amplitude(x, amplitude, peak, sigma)
        d_peak = xGaussian.der_peak(x, amplitude, d_amplitude, peak, sigma)
        d_sigma = xGaussian.der_sigma(x, amplitude, d_amplitude, peak, sigma)
        return numpy.array([d_amplitude, d_peak, d_sigma]).transpose()

    def init_parameters(self, xdata, ydata, sigma):
        """Overloaded init_parameters() method.
        """
        self.set_parameter('Amplitude', numpy.max(ydata))
        self.set_parameter('Peak', numpy.mean(xdata))
        self.set_parameter('Sigma', numpy.std(xdata))

    def fwhm(self):
        """Return the absolute FWHM of the model.
        """
        return self.SIGMA_TO_FWHM * self.sigma

    def resolution(self):
        """Return the resolution of the model, i.e., the FWHM divided by the
        peak value.
        """
        if self.peak > 0:
            return self.fwhm() / self.peak
        return 0.

    def resolution_error(self):
        """Return the error on the resolution.
        """
        if self.peak > 0 and self.parameter_error('Sigma') > 0:
            return self.resolution() * self.parameter_error('Sigma') /\
                self.parameter_value('Sigma')
        return 0.



class xFe55(xGaussian):
    """One-dimensional double gaussian model for Fe55 with 2 lines

    .. math::
      f(x; A, \\mu, \\sigma) = A e^{-\\frac{(x - \\mu)^2}{2\\sigma^2}} + cose
    """

    PARAMETER_NAMES = ('Amplitude0', 'Amplitude1', 'Peak', 'Sigma')
    PARAMETER_DEFAULT_VALUES = (1., 0.25, 0., 1.)
    PARAMETER_DEFAULT_BOUNDS = ([0., 0., -numpy.inf, 0.], [numpy.inf] * 4)
    DEFAULT_PLOTTING_RANGE = (-5., 5.)
    PEAK_SCALE = 6.490 / 5.895
    SIGMA_SCALE = numpy.sqrt(PEAK_SCALE)

    @staticmethod
    def _main_gaussian(x, amplitude0, peak, sigma):
        """
        """
        return xGaussian.value(x, amplitude0, peak, sigma)

    @staticmethod
    def _secondary_gaussian(x, amplitude1, peak, sigma):
        """
        """
        return xGaussian.value(x, amplitude1, xFe55.PEAK_SCALE * peak,
                              xFe55.SIGMA_SCALE * sigma)

    @staticmethod
    def value(x, amplitude0, amplitude1, peak, sigma):
        """Overloaded value() method.
        """
        return xGaussian.value(x, amplitude0, peak, sigma) + \
            xGaussian.value(x, amplitude1, xFe55.PEAK_SCALE * peak,
                           xFe55.SIGMA_SCALE * sigma)

    @staticmethod
    def jacobian(x, amplitude0, amplitude1, peak, sigma):
        """Overloaded jacobian() method.
        """
        peak1 = xFe55.PEAK_SCALE * peak
        sigma1 = xFe55.SIGMA_SCALE * sigma
        d_amplitude0 = xGaussian.der_amplitude(x, amplitude0, peak, sigma)
        d_amplitude1 = xGaussian.der_amplitude(x, amplitude0, peak1, sigma1)
        d_peak = xGaussian.der_peak(x, amplitude0, d_amplitude0, peak,
                                       sigma) + \
                 xGaussian.der_peak(x, amplitude1, d_amplitude1, peak1,
                                       sigma1)
        d_sigma = xGaussian.der_sigma(x, amplitude0, d_amplitude0, peak,
                                     sigma) +\
                  xGaussian.der_sigma(x, amplitude1, d_amplitude1, peak1,
                                         sigma1)
        return numpy.array([d_amplitude0, d_amplitude1,
                            d_peak, d_sigma]).transpose()

    def init_parameters(self, xdata, ydata, sigma):
        """
        """
        pass

    def plot(self, *parameters, **kwargs):
        """Overloaded plot() method.
        """
        amplitude0, amplitude1, peak, sigma = self.parameters
        xFitModelBase.plot(self, *parameters, **kwargs)
        x = numpy.linspace(self.xmin, self.xmax, 1000)
        y = self(x, *self.parameters)
        y = xGaussian.value(x, amplitude0, peak, sigma)
        plt.plot(x, y, ls='dashed')
        peak1 = xFe55.PEAK_SCALE * peak
        sigma1 = xFe55.SIGMA_SCALE * sigma
        y = xGaussian.value(x, amplitude1, peak1, sigma1)
        plt.plot(x, y, ls='dashed')


class xModulationCurveRad(xFitModelBase):

    """Modulation curve model (for fitting azimuthal distributions expressed
    in radians).

    .. math::
      f(x; A, m, \\bar\\varphi) =
      A (1 + m \\cos\\left(2(x - \\bar\\varphi)\\right))
    """

    PARAMETER_NAMES = ('Amplitude', 'Modulation', 'Phase')
    PARAMETER_DEFAULT_VALUES = (1., 1., 0.)
    PARAMETER_DEFAULT_BOUNDS = ([0., 0., -0.5 * numpy.pi],
                                [numpy.inf, 1., 0.5 * numpy.pi])
    DEFAULT_PLOTTING_RANGE = (-numpy.pi, numpy.pi)
    DEFAULT_STAT_BOX_POSITION = 'lower right'

    @staticmethod
    def value(x, amplitude, modulation, phase):
        """Overloaded value() method.
        """
        return amplitude * (1 + modulation * numpy.cos(2 * (x - phase)))

    @staticmethod
    def jacobian(x, amplitude, modulation, phase):
        """Overloaded jacobian() method.
        """
        d_amplitude = xModulationCurveRad.value(x, 1., modulation, phase)
        d_modulation = amplitude * (numpy.cos(2 * (x - phase)))
        d_phase = -2. * amplitude * modulation * numpy.sin(2 * (phase - x))
        return numpy.array([d_amplitude, d_modulation, d_phase]).transpose()

    def init_parameters(self, xdata, ydata, sigma):
        """Overloaded init_parameters() method.
        """
        self.set_parameter('Amplitude', numpy.mean(ydata))
        ymin = numpy.min(ydata)
        ymax = numpy.max(ydata)
        self.set_parameter('Modulation', (ymax - ymin) / (ymax + ymin))
        self.set_parameter('Phase', fold_angle_rad(xdata[numpy.argmax(ydata)]))


class xModulationCurveDeg(xFitModelBase):

    """Modulation curve model (for fitting azimuthal distributions expressed
    in degrees).

    .. math::
      f(x; A, m, \\bar\\varphi) =
      A (1 + m \\cos\\left(2(x - \\bar\\varphi)\\right))
    """

    PARAMETER_NAMES = ('Amplitude', 'Modulation', 'Phase')
    PARAMETER_DEFAULT_VALUES = (1., 1., 0.)
    PARAMETER_DEFAULT_BOUNDS = ([0., 0., -90.], [numpy.inf, 1., 90.])
    DEFAULT_PLOTTING_RANGE = (-180., 180.)
    DEFAULT_STAT_BOX_POSITION = 'lower right'
    DEG_TO_RAD = numpy.pi / 180.

    @staticmethod
    def value(x, amplitude, modulation, phase):
        """Overloaded value() method.

        Here we are essentially calling the modulation curve model expressed
        in radians converting to radians the input values of the independent
        variable and the phase parameter.
        """
        return xModulationCurveRad.value(numpy.radians(x), amplitude,
                                         modulation, numpy.radians(phase))

    @staticmethod
    def jacobian(x, amplitude, modulation, phase):
        """Overloaded jacobian() method.

        In addition to the degree-to-radians conversion of the inputs, here we
        have to multiply the last column of the jacobian (i.e., the phase
        derivatives) by pi / 180.
        """
        _jac = xModulationCurveRad.jacobian(numpy.radians(x), amplitude,
                                           modulation, numpy.radians(phase))
        _jac[:,2] = xModulationCurveDeg.DEG_TO_RAD * _jac[:,2]
        return _jac

    def init_parameters(self, xdata, ydata, sigma):
        """Overloaded init_parameters() method.
        """
        self.set_parameter('Amplitude', numpy.mean(ydata))
        ymin = numpy.min(ydata)
        ymax = numpy.max(ydata)
        self.set_parameter('Modulation', (ymax - ymin) / (ymax + ymin))
        self.set_parameter('Phase', fold_angle_deg(xdata[numpy.argmax(ydata)]))


class xPowerLaw(xFitModelBase):

    """Power law model.

    .. math::
      f(x; N, \\Gamma) = N x^\\Gamma
    """

    PARAMETER_NAMES = ('Normalization', 'Index')
    PARAMETER_DEFAULT_VALUES = (1., -1.)
    PARAMETER_DEFAULT_BOUNDS = ([0., -numpy.inf], [numpy.inf, numpy.inf])
    DEFAULT_PLOTTING_RANGE = (1.e-2, 1.)

    @staticmethod
    def value(x, normalization, index):
        """Overloaded value() method.
        """
        return normalization * (x**index)

    @staticmethod
    def jacobian(x, normalization, index):
        """Overloaded jacobian() method.
        """
        d_normalization = (x**index)
        d_index = numpy.log(x) * normalization * (x**index)
        return numpy.array([d_normalization, d_index]).transpose()



class xPowerLawExpCutoff(xFitModelBase):

    """Power law with an exponential cutoff.

    .. math::
      f(x; N, \\Gamma, x_0) = N x^\\Gamma e^{-\\frac{x}{x_0}}
    """

    PARAMETER_NAMES = ('Normalization', 'Index', 'Cutoff')
    PARAMETER_DEFAULT_VALUES = (1., -1., 1.)
    PARAMETER_DEFAULT_BOUNDS = ([0., -numpy.inf, 0.], [numpy.inf] * 3)
    DEFAULT_PLOTTING_RANGE = (1.e-2, 1.)

    @staticmethod
    def value(x, normalization, index, cutoff):
        """Overloaded value() method.
        """
        return normalization * (x**index) * numpy.exp(-x / cutoff)

    @staticmethod
    def jacobian(x, normalization, index, cutoff):
        """Overloaded jacobian() method.
        """
        d_normalization = (x**index) * numpy.exp(-x / cutoff)
        val = xPowerLawExpCutoff.value(x, normalization, index, cutoff)
        d_index = numpy.log(x) * val
        d_cutoff = val * x / (cutoff**2)
        return numpy.array([d_normalization, d_index, d_cutoff]).transpose()



class xPixelPha(xFitModelBase):

    """Composite model for fitting the single-pixel pulse-height distributions.

    This is essentially the sum of a gaussian centered at zero (and modeling
    the noise) and a power law with exponential cutoff (modeling the signal).

    .. math::
      f(x; A_n, \\sigma_n, N_s, \\Gamma_s, x_0) =
      A_n e^{-\\frac{x^2}{2\\sigma_n^2}} +
      N_s x^\\Gamma_s e^{-\\frac{x}{x_0}}
    """

    PARAMETER_NAMES = ('Noise amplitude', 'Noise sigma',
                       'Signal normalization', 'Signal index', 'Signal cutoff')
    PARAMETER_DEFAULT_VALUES = (1., 1., 1., -1., 1.)
    DEFAULT_PLOTTING_RANGE = (0.5, 10.)

    @staticmethod
    def _noise(x, amplitude, sigma):
        """Noise component.
        """
        return xGaussian.value(x, amplitude, 0., sigma)

    @staticmethod
    def _signal(x, normalization, index, cutoff):
        """Signal component.
        """
        return xPowerLawExpCutoff.value(x, normalization, index, cutoff)

    @staticmethod
    def value(x, noise_amplitude, noise_sigma, signal_normalization,
              signal_index, signal_cutoff):
        """Overloaded value() method.
        """
        return xPixelPha._noise(x, noise_amplitude, noise_sigma) + \
            xPixelPha._signal(x, signal_normalization, signal_index,
                             signal_cutoff)

    def plot(self, *parameters, **kwargs):
        """Overloaded plot() method.
        """
        xFitModelBase.plot(self, *parameters, **kwargs)
        x = numpy.linspace(self.xmin, self.xmax, 1000)
        y = self(x, *self.parameters)
        ymin = y.min()
        ymax = y.max()
        y = xPixelPha._noise(x, *self.parameters[:2])
        _mask = (y >= ymin) * (y <= ymax)
        plt.plot(x[_mask], y[_mask], ls='dashed')
        y = xPixelPha._signal(x, *self.parameters[2:])
        _mask = (y >= ymin) * (y <= ymax)
        plt.plot(x[_mask], y[_mask], ls='dashed')

    @staticmethod
    def jacobian(x, noise_amplitude, noise_sigma, signal_normalization,
                 signal_index, signal_cutoff):
        """Overloaded jacobian() method.
        """
        d_noise_amplitude = xGaussian.value(x, 1., 0., noise_sigma)
        d_noise_sigma = noise_amplitude * (x)**2. * d_noise_amplitude /\
                        (noise_sigma**3.)
        d_signal_normalization = (x)**signal_index *\
                                 numpy.exp(-x / signal_cutoff)
        val = xPowerLawExpCutoff.value(x, signal_normalization, signal_index,
                                          signal_cutoff)
        d_signal_index = numpy.log(x) * val
        d_signal_cutoff = val * x / (signal_cutoff**2.)
        return numpy.array([d_noise_amplitude, d_noise_sigma,
                            d_signal_normalization, d_signal_index,
                            d_signal_cutoff]).transpose()



class xExponential(xFitModelBase):

    """Exponential model.

    .. math::
      f(x; N, \\alpha) = N e^{\\alpha x}
    """

    PARAMETER_NAMES = ('Normalization', 'Index')
    PARAMETER_DEFAULT_VALUES = (1., -1.)
    PARAMETER_DEFAULT_BOUNDS = ([0., -numpy.inf], [numpy.inf]*2)
    DEFAULT_PLOTTING_RANGE = (0., 1.)

    @staticmethod
    def value(x, normalization, exponent):
        """Overloaded value() method.
        """
        return normalization * numpy.exp(exponent * x)

    @staticmethod
    def jacobian(x, normalization, exponent):
        """Overloaded jacobian() method.
        """

        d_normalization = numpy.exp(exponent * x)
        d_exponent = normalization * x * numpy.exp(exponent * x)
        return numpy.array([d_normalization, d_exponent]).transpose()



class xExponentialOffset(xFitModelBase):

    """Exponential model.

    .. math::
      f(x; A, N, \\alpha) = A + N e^{\\alpha x}
    """

    PARAMETER_NAMES = ('Offset', 'Normalization', 'Index')
    PARAMETER_DEFAULT_VALUES = (0., 1., -1.)
    PARAMETER_DEFAULT_BOUNDS = ([-numpy.inf, 0., -numpy.inf], [numpy.inf]*3)
    DEFAULT_PLOTTING_RANGE = (0., 1.)

    @staticmethod
    def value(x, offset, normalization, exponent):
        """Overloaded value() method.
        """
        return offset + normalization * numpy.exp(exponent * x)

    @staticmethod
    def jacobian(x, offset, normalization, exponent):
        """Overloaded jacobian() method.
        """
        d_offset = numpy.full((len(x),), 1.)
        d_normalization = numpy.exp(exponent * x)
        d_exponent = normalization * x * numpy.exp(exponent * x)
        return numpy.array([d_offset, d_normalization, d_exponent]).transpose()



class xLogNormal(xFitModelBase):

    """Log-normal fitting model.

    See
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html
    for more details on the implementation.

    Note that our parametrization is done in terms of the mean and sigma of the
    distribution, in addition to the shape parameter.
    """

    PARAMETER_NAMES = ('Normalization', 'Mean', 'Sigma', 'Shape')
    PARAMETER_DEFAULT_VALUES = (1., 0., 1., 0.02)
    DEFAULT_PLOTTING_RANGE = (-5., 5.)

    @staticmethod
    def value(x, norm, mu, sigma, shape):
        """Overloaded value() method.
        """
        p = numpy.exp(shape * shape)
        scale = sigma / numpy.sqrt(p * (p - 1.))
        loc = mu - scale * numpy.sqrt(p)
        return norm * scipy.stats.lognorm.pdf(x, shape, loc, scale)


class xGeneralizedGaussian(xFitModelBase):

    """Generalized gaussian fitting model.

    See
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gennorm.html
    for implementation details.
    """

    PARAMETER_NAMES = ('Normalization', 'Mean', 'Sigma', 'Beta')
    PARAMETER_DEFAULT_VALUES = (1., 0., 1., 1.)
    DEFAULT_PLOTTING_RANGE = (-5., 5.)

    @staticmethod
    def value(x, norm, mean, sigma, beta):
        """Overloaded value() method.
        """
        return norm * scipy.stats.gennorm.pdf(x, beta, mean, sigma)



class xHat(xFitModelBase):

    """Fit model representing a flat distribution truncated with an error
    function at both ends.

    The ASYMMETRY class member controls the slope of the central part of the model.
    """

    PARAMETER_NAMES = ('Normalization', 'X1', 'Sigma1', 'X2', 'Sigma2')
    PARAMETER_DEFAULT_VALUES = (1., 0.2, 0.05, 0.8, 0.05)
    DEFAULT_PLOTTING_RANGE = (0., 1.)
    ASYMMETRY = 0.

    @staticmethod
    def value(x, norm, x1, sigma1, x2, sigma2):
        """Overloaded value() method.
        """
        z1 = (x - x1) / sigma1
        z2 = (x - x2) / sigma2
        dx = (x2 - x1)
        val = norm * 0.25 * (1. + scipy.special.erf(z1)) * (1. + scipy.special.erf(-z2)) / dx
        if xHat.ASYMMETRY != 0:
            val *= 1. - xHat.ASYMMETRY * (x - x1) / dx
        return val


class xLorentzian(xFitModelBase):

    """Lorentzian function
    """

    PARAMETER_NAMES = ('Normalization', 'Peak', 'HWHM' )
    PARAMETER_DEFAULT_VALUES = (1., 0., 1.)
    PARAMETER_DEFAULT_BOUNDS = ([0., -numpy.inf, 0.], [numpy.inf, numpy.inf, numpy.inf])
    DEFAULT_PLOTTING_RANGE = (-1., 1.)

    @staticmethod
    def value(x, normalization, peak, hwhm):
        """Overloaded value() method.
        """
        return normalization * hwhm / ((x - peak)**2. + hwhm**2.)

    @staticmethod
    def jacobian(x, normalization, peak, hwhm):
        """Overloaded jacobian() method.
        """
        d_x = x - peak
        denom = d_x**2. + hwhm**2.
        d_normalization = hwhm / denom
        d_peak = 2. * normalization * hwhm * d_x / denom**2.
        d_hwhm = normalization * (d_x**2.  - hwhm**2.)/ denom**2.
        return numpy.array([d_normalization, d_peak, d_hwhm]).transpose()

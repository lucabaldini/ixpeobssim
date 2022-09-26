#!/usr/bin/env python
#
# Copyright (C) 2018--2021, the ixpeobssim team.
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
from scipy.optimize import curve_fit

from ixpeobssim.utils.environment import SCIPY_VERSION
import ixpeobssim.core.modeling
from ixpeobssim.utils.logging_ import logger

# pylint: disable=invalid-name


USE_ABSOLUTE_SIGMA = True


def fit(model, xdata, ydata, p0=None, sigma=None, xmin=-numpy.inf,
        xmax=numpy.inf, absolute_sigma=USE_ABSOLUTE_SIGMA, check_finite=True,
        method=None, verbose=True, **kwargs):
    """Lightweight wrapper over the ``scipy.optimize.curve_fit()`` function
    to take advantage of the modeling facilities. More specifically, in addition
    to performing the actual fit, we update all the model parameters so that,
    after the fact, we do have a complete picture of the fit outcome.

    Parameters
    ----------
    model : :py:class:`ixpeobssim.core.modeling.FitModelBase` instance callable
        The model function, f(x, ...).  It must take the independent
        variable as the first argument and the parameters to fit as
        separate remaining arguments.

    xdata : array_like
        The independent variable where the data is measured.

    ydata : array_like
        The dependent data --- nominally f(xdata, ...)

    p0 : None, scalar, or sequence, optional
        Initial guess for the parameters. If None, then the initial
        values will all be 1.

    sigma : None or array_like, optional
        Uncertainties in `ydata`. If None, all the uncertainties are set to
        1 and the fit becomes effectively unweighted.

    xmin : float
        The minimum value for the input x-values.

    xmax : float
        The maximum value for the input x-values.

    absolute_sigma : bool, optional
        If True, `sigma` is used in an absolute sense and the estimated
        parameter covariance `pcov` reflects these absolute values.
        If False, only the relative magnitudes of the `sigma` values matter.
        The returned parameter covariance matrix `pcov` is based on scaling
        `sigma` by a constant factor. This constant is set by demanding that the
        reduced `chisq` for the optimal parameters `popt` when using the
        *scaled* `sigma` equals unity.

    method : {'lm', 'trf', 'dogbox'}, optional
        Method to use for optimization.  See `least_squares` for more details.
        Default is 'lm' for unconstrained problems and 'trf' if `bounds` are
        provided. The method 'lm' won't work when the number of observations
        is less than the number of variables, use 'trf' or 'dogbox' in this
        case.

    verbose : bool
        Print the model if True.

    kwargs
        Keyword arguments passed to `leastsq` for ``method='lm'`` or
        `least_squares` otherwise.
    """
    # Select data based on the x-axis range passed as an argument.
    _mask = numpy.logical_and(xdata >= xmin, xdata <= xmax)
    xdata = xdata[_mask]
    ydata = ydata[_mask]
    if len(xdata) <= len(model.parameters):
        raise RuntimeError('Not enough data to fit (%d points)' % len(xdata))
    if isinstance(sigma, numpy.ndarray):
        sigma = sigma[_mask]
    # If the model has a Jacobian defined, go ahead and use it.
    try:
        jac = model.jacobian
    except:
        jac = None
    # If we are not passing default starting points for the model parameters,
    # try and do something sensible.
    if p0 is None:
        model.init_parameters(xdata, ydata, sigma)
        p0 = model.parameters
        if verbose:
            logger.debug('%s parameters initialized to %s.', model.name(), p0)
    # If sigma is None, assume all the errors are 1. (If we don't do this,
    # the code will crash when calculating the chisquare.
    if sigma is None:
        sigma = numpy.full((len(ydata), ), 1.)
    # Do the actual fit---note the horrible hack to support scipy versions prior
    # to 0.18, not supporting a user-defined jacobian.
    if SCIPY_VERSION < '0.18':
        logger.warning('Scipy version < 0.18 found, ignoring Jacobian...')
        popt, pcov = curve_fit(model, xdata, ydata, p0, sigma, absolute_sigma,
                               check_finite, model.bounds, method, **kwargs)
    else:
        popt, pcov = curve_fit(model, xdata, ydata, p0, sigma, absolute_sigma,
                               check_finite, model.bounds, method, jac,
                               **kwargs)
    # Update the model parameters.
    model.set_plotting_range(xdata.min(), xdata.max())
    model.parameters = popt
    model.covariance_matrix = pcov
    model.chisq = (((ydata - model(xdata))/sigma)**2).sum()
    model.ndof = len(ydata) - len(model)
    if verbose:
        print(model)
    return model


def fit_histogram(model, histogram, p0=None, sigma=None, xmin=-numpy.inf,
                  xmax=numpy.inf, absolute_sigma=USE_ABSOLUTE_SIGMA,
                  check_finite=True, method=None, verbose=True, **kwargs):
    """Fit a histogram to a given model.

    This is basically calling :meth:`ixpeobssim.core.fitting.fit` with some
    pre-processing to turn the histogram bin edges and content into
    x-y data. Particularly, the bin centers are taken as the independent
    data series, the bin contents are taken as the dependent data saries,
    and the square root of the counts as the Poisson error.

    For additional parameters look at the documentation of the
    :meth:`ixpeobssim.core.fitting.fit`

    Parameters
    ----------
    model : :py:class:`ixpeobssim.core.modeling.FitModelBase` instance or
        callable
        The fit model.

    histogram : ixpeHistogram1d instance
        The histogram to be fitted.

    Warning
    -------
    We're not quite doing the right thing, here, as we should integrate
    the model within each histogram bin and compare that to the counts,
    but this is not an unreasonable first-order approximation. We might want
    to revise this, especially since we can probably provide an analytic
    integral for most of the model we need.
    """
    assert histogram.num_axes == 1
    _mask = (histogram.content > 0)
    xdata = histogram.bin_centers(0)[_mask]
    ydata = histogram.content[_mask]
    if sigma is None:
        sigma = numpy.sqrt(ydata)
    return fit(model, xdata, ydata, p0, sigma, xmin, xmax, absolute_sigma,
               check_finite, method, verbose, **kwargs)


def fit_modulation_curve(histogram, p0=None, sigma=None, xmin=-numpy.inf,
                         xmax=numpy.inf, absolute_sigma=USE_ABSOLUTE_SIGMA,
                         check_finite=True, method=None, verbose=True,
                         degrees=True, **kwargs):
    """Fit a modulation curve.

    For all the other parameters look at the documentation of the
    :meth:`ixpeobssim.core.fitting.fit_histogram`
    """
    if degrees:
        model = ixpeobssim.core.modeling.xModulationCurveDeg()
    else:
        model = ixpeobssim.core.modeling.xModulationCurveRad()
    return fit_histogram(model, histogram, p0, sigma, xmin, xmax,
                         absolute_sigma, check_finite, method, verbose,
                         **kwargs)


def fit_gaussian_iterative(histogram, p0=None, sigma=None, xmin=-numpy.inf,
                           xmax=numpy.inf, absolute_sigma=USE_ABSOLUTE_SIGMA,
                           check_finite=True, method=None, verbose=True,
                           num_sigma_left=2., num_sigma_right=2.,
                           num_iterations=2, **kwargs):
    """Fit the core of a gaussian histogram within a given number of sigma
    around the peak.

    This function performs a first round of fit to the data and then
    repeats the fit iteratively limiting the fit range to a specified
    interval defined in terms of deviations (in sigma) around the peak.

    For additional parameters look at the documentation of the
    :meth:`ixpeobssim.core.fitting.fit_histogram`

    Parameters
    ----------
    num_sigma_left : float
        The number of sigma on the left of the peak to be used to define the
        fitting range.

    num_sigma_right : float
        The number of sigma on the right of the peak to be used to define the
        fitting range.

    num_iterations : int
        The number of iterations of the fit.
    """
    model = ixpeobssim.core.modeling.xGaussian()
    fit_histogram(model, histogram, p0, sigma, xmin, xmax, absolute_sigma,
                  check_finite, method, verbose, **kwargs)
    for i in range(num_iterations):
        xmin = model.peak - num_sigma_left * model.sigma
        xmax = model.peak + num_sigma_right * model.sigma
        try:
            fit_histogram(model, histogram, p0, sigma, xmin, xmax,
                          absolute_sigma, check_finite, method, verbose,
                          **kwargs)
        except RuntimeError as e:
            raise RuntimeError('%s after %d iteration(s)' % (e, i + 1))
    return model



def _analytical_fit_base(x, y):
    """Basic common routine for performing analytical fits.

    The basic idea is to have a framework to perform simple, un-weighted least-squares
    fits (e.g., with linear or power-law models) that have a closed-form solution
    and can be performed analytically. Note this is not computing fit errors, but
    on the plus side it offers limited support for vectorized fits, which is
    sometimes useful to perform extrapolation of noisy data. More precisely,
    we support three cases:

    * x and y are two one dimensional arrays of the same shape, in which case
      we provide scalar fit parameters;
    * x is a 1-dimensional array with shape (N,), and y is a 2-dimensional
      arrays with shape (N, M), in which case we interpret the situation as M
      fits where the independent data array is the same for all the data arrays
      of the dependent variable, and an array of values of lenght M is returned
      for each parameter;
    * x and y are both 2-dimensional arrays with shape (N, M), which is
      similar to the previous case, except for the fact that the M data
      arrays for the independent variable are all different.

    Arguments
    ---------
    x : array_like
        The x data, can be either one- (N) or two- (N, M) dimensional.

    y : array_like
        The y data, can be either one- (N) or two- (N, M) dimensional.
    """
    # Make sure the input data are in the form on numpy arrays.
    if not isinstance(x, numpy.ndarray):
        x = numpy.array(x)
    if not isinstance(y, numpy.ndarray):
        y = numpy.array(y)
    # Handle the common case where we pass a single array for x and multiple arrays
    # for the y---in this case we simply tile the x with enough identical copies
    # to match the y shape.
    if x.shape != y.shape:
        x = numpy.tile(x, (y.shape[1], 1)).T
    # And, at this point, the shapes of x and y should really match.
    assert x.shape == y.shape
    n = x.shape[0]
    return x, y, n


def linear_analytical_fit(x, y, dy=None):
    """Simple vectorized, analytical linear fit.

    See https://mathworld.wolfram.com/LeastSquaresFitting.html
    """
    x, y, _ = _analytical_fit_base(x, y)
    if dy is None:
        dy = numpy.full(y.shape, 1.)
    w = 1. / dy**2.
    s0 = w.sum(axis=0)
    sx = (x * w).sum(axis=0)
    sx2 = (x**2. * w).sum(axis=0)
    sy = (y * w).sum(axis=0)
    sxy = (x * y * w).sum(axis=0)
    D = s0 * sx2 - sx**2.
    m = (s0 * sxy - sx * sy) / D
    q = (sy * sx2 - sx * sxy) / D
    return m, q


def power_law_analytical_fit(x, y):
    """Simple vectorized, analytical linear fit.

    See https://mathworld.wolfram.com/LeastSquaresFitting.html
    """
    x, y, n = _analytical_fit_base(x, y)
    logx, logy = numpy.log(x), numpy.log(y)
    sx = logx.sum(axis=0)
    sx2 = (logx**2.).sum(axis=0)
    sy = logy.sum(axis=0)
    sxy = (logx * logy).sum(axis=0)
    index = (n * sxy - sx * sy) / (n * sx2 - sx**2.)
    norm = numpy.exp((sy - index * sx) / n)
    return norm, index

#!/usr/bin/env python
#
# Copyright (C) 2015--2019, the ixpeobssim team.
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

"""Spline utility module, building on top of the scipy.interpolate modules.
"""

from __future__ import print_function, division

from numbers import Number

import matplotlib
import numpy
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import RectBivariateSpline

from ixpeobssim.utils.matplotlib_ import plt, last_line_color, setup_gca
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.misc import pairwise, pairwise_enum

# pylint: disable=invalid-name, too-many-arguments

def interpolate(xa, ya, xb, yb, x):
    """Simple two-point linear interpolation/extrapolation.

    This is a convenience function that is pretty much only used in the
    optimize_grid_linear() below.
    """
    return ya + (yb - ya) / (xb - xa) * (x - xa)

def optimize_grid_linear(x, y, tolerance=1e-4):
    """Optimize a pair of (x, y) arrays for the corresponding spline
    definition.

    This loops over the input arrays and removes unnecessary data points
    to minimize the length of the arrays necessary to the spline definition.

    Args
    ----
    x : array
        The input x-array.

    y : array
        The input y-array.

    tolerance : float
        The maximum relative difference between the generic yi value and the\
        estrapolation of the two previous optimized data points for the point\
        i to be removed.
    """
    assert len(x) == len(y)
    logger.info('Optimizing grid with %d starting points...', len(x))
    # Start a new series with the first two points of the input arrays.
    _x = [x[0], x[1]]
    _y = [y[0], y[1]]
    # Loop over the points 3 ... (N - 1).
    for (_xi, _yi) in list(zip(x, y))[2:-1]:
        # Extrapolate the last two points of the new series to xi and
        # see how far we are from the actual yi.
        delta = interpolate(_x[-2], _y[-2], _x[-1], _y[-1], _xi) - _yi
        if abs(delta / _yi) > tolerance:
            # If the difference is larger than the tolerance, add a point.
            # (This has the drawback that we tend to add pairs of point at
            # each change of slope.)
            _x.append(_xi)
            _y.append(_yi)
            # Interpolate the points last and (last - 2) to (last - 1).
            delta = interpolate(_x[-3], _y[-3], _x[-1], _y[-1], _x[-2]) - _y[-2]
            if abs(delta / _y[-2]) < tolerance:
                # If the penultimate point was not necessary, remove it.
                _x.remove(_x[-2])
                _y.remove(_y[-2])
    # Append the last point of the original array to the list.
    _x.append(x[-1])
    _y.append(y[-1])
    _x, _y = numpy.array(_x), numpy.array(_y)
    logger.info('Done, %d points remaining.', len(_x))
    return _x, _y


class xUnivariateSpline(UnivariateSpline):

    """Wrapper around scipy's UnivariateSplineClass.

    The basic idea is to keep track of the original arrays passed to the
    interpolator and to support arithmetic operations and plotting. We also
    allow the user to supply optional arguments to control the ranges and
    specify labels (e.g., names and units) for the quantities involved.
    We essentially maintain the original signature of the function, adding
    two optional arguments (xlabel and ylabel) at the end.

    Note
    ----
    Note that the interface to the base class has changed from numpy 0.14.
    An `ext` argument can be passed to the constructor starting with scipy
    0.15 to control the extrapolation behavior and a `check_finite` argument is
    available in 0.16 to avoid `nans` in the input data. We currently do not
    expose either one and rely on the default parameters for backward
    compatibility.

    Args
    ----
    x : (N,) array_like
        Input x values for the spline.

    y : (N,) array_like
        Input y values for the spline.

    w : (N,) array_like, optional
        Weights for spline fitting.

    bbox : (2,) array_like, optional
        Bounduary of the approximation interval. .

    k : int
        The order of the spline (<=5, default 3, i.e., a cubic spline).

    s : float or None
        The spline smoothing factor (o means no smoothing).

    ext : int or string, optional
        Controls the extrapolation mode for elements not in the interval defined
        by the knot sequence.

    xlabel: str, optional
        A text label for the quantity on the x-axis.

    ylabel: str, optional
        A text label for the quantity on the y-axis.

    log : bool, default False
        Flag to initialize the underline scipy spline in logarithmic space.
    """

    def __init__(self, x, y, w=None, bbox=None, k=3, s=None, ext=0,
                 xlabel=None, ylabel=None, log=False):
        """Constructor.
        """
        # Make sure the input vectors have the same lengths.
        assert len(x) == len(y)
        # numpy.unique() is returning a sorted copy of the unique values of the
        # x arrays, so this is effectively sorting x.
        self.x, _index = numpy.unique(x, return_index=True)
        # If some of the values were not unique, give up.
        assert len(self.x) == len(x)
        # Make sure the input array is long enough for a spline of order k.
        if len(self.x) <= k:
            abort('x spline array (%s) too short for k = %d' % (self.x, k))
        # Ok, now that the x values are fine we need to take the y elements
        # in the right order.
        self.y = y[_index]
        # This is done so that we don't have to pass mutable argument defaults.
        if bbox is None:
            bbox = [None, None]
        # Bookkeeping for the labels.
        self.xlabel = xlabel
        self.ylabel = ylabel
        # And, finally: initialize the actual scipy spline.
        if not log:
            args = (self.x, self.y, w, bbox, k, s, ext)
        else:
            args = (numpy.log10(self.x), numpy.log10(self.y), w, bbox, k, s, ext)
        UnivariateSpline.__init__(self, *args)

    def xmin(self):
        """Return the minimum of the underlying x-array.

        Note this works because the input array is sorted.
        """
        return self.x[0]

    def xmax(self):
        """Return the maximum of the underlying x-array.

        Note this works because the input array is sorted.
        """
        return self.x[-1]

    def _xunion(self, other):
        """Small convenience function to merge the x values of two splines.

        This is used in the overloading of the comparison operators.
        """
        return numpy.union1d(self.x, other.x)

    def __mul__(self, other):
        """Overloaded multiplication operator.

        Note we do not attempt to calculate a sensible value for the y-label,
        in this case. The user can adjust it after the fact.
        """
        x = self._xunion(other)
        y = self(x) * other(x)
        return self.__class__(x, y, xlabel=self.xlabel)

    def __truediv__(self, other):
        """Overloaded division operator.

        Note we do not attempt to calculate a sensible value for the y-label,
        in this case. The user can adjust it after the fact.
        """
        x = self._xunion(other)
        x = x[other(x) != 0]
        y = self(x) / other(x)
        return self.__class__(x, y, xlabel=self.xlabel)

    def __add__(self, other):
        """Overloaded sum operator.
        """
        x = self._xunion(other)
        y = self(x) + other(x)
        return self.__class__(x, y, xlabel=self.xlabel, ylabel=self.ylabel)

    def __sub__(self, other):
        """Overloaded sum operator.
        """
        x = self._xunion(other)
        y = self(x) - other(x)
        return self.__class__(x, y, xlabel=self.xlabel, ylabel=self.ylabel)

    def __len__(self):
        """Return the lenght of the arrays used to construct the spline.
        """
        return len(self.x)

    def scale(self, scale_factor, **kwargs):
        """Return a different spline whose y values are scaled by a given factor
        wrt the original one.
        """
        x = numpy.copy(self.x)
        y = numpy.copy(self.y) * scale_factor
        kwargs.setdefault('xlabel', self.xlabel)
        kwargs.setdefault('ylabel', self.ylabel)
        return self.__class__(x, y, **kwargs)

    def inverse(self, xmin=None, xmax=None, **kwargs):
        """Calculate the inverse of the spline.

        Note that the xmin and xmax arguments allow to invert the spline in
        a subset of its values.

        Arguments
        ---------
        xmin : float
            The minimum value for the x-axis of the inverse spline (or y-axis of
            the original one).

        xmax : float
            The maximum value for the x-axis of the inverse spline (or y-axis of
            the original one).
        """
        # Swap the x and y arrays, making sure to get rid of the duplicated
        # y values, if any.
        x, mask = numpy.unique(self.y, return_index=True)
        y = self.x[mask]
        # Now restrict to the target interval, if necessary.
        if xmin is None:
            xmin = self.y.min()
        if xmax is None:
            xmax = self.y.max()
        mask = numpy.logical_and(x >= xmin, x <= xmax)
        x = x[mask]
        y = y[mask]
        kwargs.setdefault('xlabel', self.ylabel)
        kwargs.setdefault('ylabel', self.xlabel)
        return self.__class__(x, y, **kwargs)

    def norm(self):
        """Return the integral over the entire spline domain.
        """
        return self.integral(self.xmin(), self.xmax())

    def build_cdf(self, ylabel='cdf', spline_class=None):
        """Create the cumulative distribution function.

        Args
        ----
        spline_class : class
            The specific spline class (e.g., xInterpolatedUnivariateSpline)
            used to build the cdf

        ylabel : str
            The (optional) label for the y axis

        Note
        ----
        At this point the dafault behavior is to use a linear spline for the
        cdf by default. This is mainly for historical reasons, but maybe
        a better option would be to use a cubic spline, or to delegate
        this to call self.__class__?

        Warning
        -------
        Re-enable the check on y > 0!
        """
        # For the cdf to make sense, the input pdf must be positive-defined.
        #assert (self.y >= 0.).all()
        if spline_class is None:
            spline_class = xInterpolatedUnivariateSplineLinear
        _xmin = self.xmin()
        # Loop over the underlying x array and calculate the partial integrals.
        _y = [self.integral(_xmin, _xbar) for _xbar in self.x]
        # Tranform the list into a numpy array and divide by the last
        # element to that cdf(xmax) = 1.
        _y = numpy.array(_y) / _y[-1]
        # Build the actual spline object to be returned.
        return spline_class(self.x, _y, xlabel=self.xlabel, ylabel=ylabel)

    def build_ppf(self, xlabel='q', spline_class=None):
        """Create the percent point function (or inverse of the cdf).

        See the docstring for the build_cdf() method for more details.
        """
        # For the ppf to make sense, the input pdf must be positive-defined.
        #assert (self.y >= 0.).all()
        if spline_class is None:
            spline_class = xInterpolatedUnivariateSplineLinear
        _xmin = self.xmin()
        # Loop over the underlying x array and calculate the partial integrals.
        # Note this time this is assigned to the x axis!
        _x = numpy.array([self.integral(_xmin, _xbar) for _xbar in self.x])
        # Need an explicit extra loop to remove the non-unique values from the
        # resulting array.
        _x, _mask = numpy.unique(_x, return_index=True)
        # Stretch the x-axis of the output spline to the [0, 1] interval.
        _x /= _x[-1]
        _y = self.x[_mask]
        # Build the actual spline object to be returned.
        return spline_class(_x, _y, xlabel=xlabel, ylabel=self.xlabel)

    def plot(self, num_points=1000, overlay=False, logx=False, logy=False,
             scale=1., offset=0., grids=False, **kwargs):
        """Plot the spline.

        Args
        ----
        num_points : int, optional
            The number of sampling points to be used to draw the spline.

        overlay : bool, optional
            If True, the original arrays passed to the spline are overlaid.

        logx : bool, optional
            If True, the spline is sampled and plotted with the log scale on the
            x axis.

        logy : bool, optional
            If True, the spline is plotted with the log scale on the y axis.

        scale : float, optional
            Optional scale factor for plotting (useful for overlaying splines
            with different ranges on the y axis).

        offset : float, optional
            Optional offset for plotting (useful for overlaying splines with
            different ranges on the y axis).

        **kwargs : dict
            Keyword arguments passed directly to the matplotlib.plot() method.
        """
        # Choose a proper grid on the x axis for plotting.
        if logx:
            _x = numpy.logspace(numpy.log10(self.xmin()),
                                numpy.log10(self.xmax()), num_points)
        else:
            _x = numpy.linspace(self.xmin(), self.xmax(), num_points)
        # Evaluate the spline on the input grid, applying the scale factor and
        # the offset, if needed.
        _y = scale * self(_x) + offset
        # Plot the thing.
        plt.plot(_x, _y, '-', **kwargs)
        if overlay:
            kwargs['color'] = last_line_color()
            kwargs['label'] = None
            plt.plot(self.x, self.y, 'o', **kwargs)
        # Set up the axes labels and scales.
        setup_gca(self.xlabel, self.ylabel, logx=logx, logy=logy, grids=grids)



class xInterpolatedUnivariateSpline(xUnivariateSpline):

    """Interpolated spline (i.e., passing through all the input points).
    """

    def __init__(self, x, y, w=None, bbox=None, k=3, ext=0, xlabel=None, ylabel=None):
        """Constructor.
        """
        xUnivariateSpline.__init__(self, x, y, w, bbox, k, 0, ext, xlabel, ylabel)



class xInterpolatedPiecewiseUnivariateSpline(xInterpolatedUnivariateSpline):

    """Pieacewise interpolated spline.

    In addition to all the regular xInterpolatedUnivariateSpline parameters,
    this takes an additional breakpoints argument that allows to split the
    spline creation and evaluation into independent pieces. This is achieved
    by stories a list of independent splines under the hood.
    """

    def __init__(self, x, y, breakpoints, w=None, bbox=None, k=3, ext=0, xlabel=None, ylabel=None):
        """Constructor.
        """
        if isinstance(breakpoints, Number):
            breakpoints = (breakpoints, )
        self.breakpoints = numpy.array(breakpoints, dtype=float)
        args = w, bbox, k, ext
        self.x = x
        self.y = y
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.__splines = []
        for min_, max_ in pairwise(self.__split_points(self.x)):
            _x = x[min_:max_]
            _y = y[min_:max_]
            self.__splines.append(xInterpolatedUnivariateSpline(_x, _y, *args))

    def __split_points(self, x):
        """Return the split points for a given array of x.
        """
        splits = numpy.searchsorted(x, self.breakpoints)
        try:
            # Multiple splits?
            splits = tuple(splits)
        except TypeError:
            # Catch the case where we have a single split, and numpy.searchsorted()
            # returns a scalar.
            splits = (splits, )
        return (0, ) + splits + (len(x), )

    def __call__(self, x):
        """Overloaded method.

        The slight complication, here, is due to the fact that we're trying to
        support both scalars and numpy arrays.
        """
        # If x is a scalar,we have to bisect the breakpoints with x and find the
        # index of the appropriate spline.
        if isinstance(x, Number):
            return self.__splines[numpy.searchsorted(self.breakpoints, x)](x)
        # Otherwise we have to do the opposite---split the input x array in the
        # proper pieces and evaluate the splines separately, concatenating the
        # results.
        vals = numpy.array([], dtype=float)
        for i, (min_, max_) in pairwise_enum(self.__split_points(x)):
            vals = numpy.append(vals, self.__splines[i](x[min_:max_]))
        return vals



class xInterpolatedUnivariateSplineLinear(xUnivariateSpline):

    """Interpolate linear (k=1) spline.
    """

    def __init__(self, x, y, w=None, bbox=None, ext=0, xlabel=None, ylabel=None):
        """ Constructor.
        """
        xUnivariateSpline.__init__(self, x, y, w, bbox, 1, 0, ext, xlabel, ylabel)



class xInterpolatedUnivariateLogSpline(xUnivariateSpline):

    """Poor man's attempt at a spline in logarithmic space.

    The class inherits from both xUnivariateSplineBase and scipy's
    InterpolatedUnivariateSpline, and the basic strategy, here, is twofold:
    we initialize xUnivariateSplineBase from the physical x and y values,
    but we construct an actual interpolator in logarithmic space under the
    hood, by passing log10(x) and log10(y) to the InterpolatedUnivariateSpline
    constructor. We then overload the __call__() method by tranforming the
    argument into its logarithm, performing the interpolation in logarithmic
    space, and raising the result to the power 10.

    The class comes with some pretty obvious shortcomings, among which the
    difficulty of implementing sensible derivatives and integrals. As far as
    the latter is concerned, we proceed by brute force and create a second
    spline, this time in linear space, with a relatively large number of points
    to guarantee the accuracy of the integral. Admittedly, this is not very
    elegant nor very roboust. For the derivative issue, we simply do not
    provide the nu argument in the overloaded __call__ method. For all the
    other scipy's spline method that is not implemented in this context, we
    just overload it and raise a NotImplementedError exception.

    At some point we even considered removing this class, but the truth is
    it is handy when one needs a quick way to extrapolate in log-log space.
    So the class is remaning, with all the caveates above.
    """

    def __init__(self, x, y, w=None, bbox=None, k=3, ext=0, xlabel=None, ylabel=None):
        """Constructor.
        """
        args = (x, y, w, bbox, k, 0, ext, xlabel, ylabel)
        xUnivariateSpline.__init__(self, *args, log=True)
        self.__linear_spline = None

    # pylint: disable=signature-differs
    def __call__(self, x):
        """Overloaded call method.

        Note the signature is different from that of the parent, i.e.,
        InterpolatedUnivariateSpline.__call__(x, nu=0, ext=None)
        because the do not support derivatives, at this point.
        """
        return numpy.power(10., UnivariateSpline.__call__(self, numpy.log10(x)))

    def antiderivative(self, n=1):
        """Overloaded method.
        """
        raise NotImplementedError

    def derivative(self, n=1):
        """Overloaded method.
        """
        raise NotImplementedError

    def derivatives(self, x):
        """Overloaded method.
        """
        raise NotImplementedError

    def roots(self):
        """Overloaded method.
        """
        raise NotImplementedError

    def __build_integral_spline(self, num_points=1000, k=3):
        """Build the underlying univariate linear interpolated spline to be
        used to evaluate integrals.

        Warning
        -------
        We are not exposing the num_points and k parameters anywhere in the
        methods calling this one, so they are effectively fixed. If necessary
        we might need to revise this.
        """
        _logxmin = numpy.log10(self.xmin())
        _logxmax = numpy.log10(self.xmax())
        _x = numpy.logspace(_logxmin, _logxmax, num_points)
        _y = self(_x)
        fmt = dict(xlabel=self.xlabel, ylabel=self.ylabel)
        return xInterpolatedUnivariateSpline(_x, _y, k=k, **fmt)

    def integral(self, a, b):
        """Overloaded integral method.

        The integral spline is calculated and cached the first time this method
        is called.
        """
        if self.__linear_spline is None:
            self.__linear_spline = self.__build_integral_spline()
        return self.__linear_spline.integral(a, b)



class xInterpolatedUnivariateLogSplineLinear(xInterpolatedUnivariateLogSpline):

    """Subclass of xInterpolatedUnivariateLogSpline with k=1.
    """

    def __init__(self, x, y, w=None, bbox=None, ext=0, xlabel=None, ylabel=None):
        """Constructor.
        """
        args = (x, y, w, bbox, 1, ext, xlabel, ylabel)
        xInterpolatedUnivariateLogSpline.__init__(self, *args)

    def antiderivative(self, n=1):
        """Overloaded method.
        """
        raise NotImplementedError

    def derivative(self, n=1):
        """Overloaded method.
        """
        raise NotImplementedError

    def derivatives(self, x):
        """Overloaded method.
        """
        raise NotImplementedError

    def roots(self):
        """Overloaded method.
        """
        raise NotImplementedError



class xInterpolatedBivariateSpline(RectBivariateSpline):

    """Bivariate interpolated spline on a rectangular grid.

    This is somewhat similar in spirit to the corresponding univariate base
    class, except that the additional functionalities are, for the moment,
    limited to book-keeping and plotting facilities.

    One handy facility that this class provides is that of allowing both
    an array like of shape (x.size, y.size) and a function as the z argument.
    This makes it possible to construct slides by passing either (x, y, z), as
    the signature of the underlying scipy's interpolator, and (x, y, f(x, y)),
    which is easier and more intuitive in most cases.

    Args
    ----
    x : array_like
        Input x values (assumed to be sorted).

    y : array_like
        Input y values (assumed to be sorted).

    z : array_like of shape (x.size, y.size) or callable
        Input z values, with shape (x.size,y.size).

    bbox : array_like, optional
        The boundary of the approximation domain.

    kx, ky : int
        The spline order on the two axes.

    xlabel: str, optional
        A text label for the quantity on the x-axis.

    ylabel: str, optional
        A text label for the quantity on the y-axis.

    zlabel: str, optional
        A text label for the quantity on the z-axis.
    """

    def __init__(self, x, y, z, bbox=None, kx=3, ky=3, xlabel=None, ylabel=None,
                 zlabel=None):
        """Constructor.
        """
        self.x = x
        self.y = y
        # If the z argument is callable, evaluate it on the appropriate grid
        # to build the spline (this makes it easy to construct interpolators
        # by passing x, y and f(x, y)).
        if hasattr(z, '__call__'):
            z = z(*self.transposed_meshgrid(x, y))
        self.z = z
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel
        if bbox is None:
            bbox = [None, None, None, None]
        RectBivariateSpline.__init__(self, x, y, z, bbox, kx, ky, s=0)

    @staticmethod
    def transposed_meshgrid(x, y):
        """Return a suitable meshgrid to evaluate a function in the spline
        reference system.

        This is needed because apparently meshgrid is using the matlab
        convention for the indices, which is inconsistent with most of the
        other interfaces in scipy. See:
        https://github.com/scipy/scipy/issues/3164
        """
        _x, _y = numpy.meshgrid(x, y)
        return _x.T, _y.T

    def xmin(self):
        """Return the minimum of the underlying x-array.
        """
        return self.x[0]

    def xmax(self):
        """Return the maximum of the underlying x-array.
        """
        return self.x[-1]

    def ymin(self):
        """Return the minimum of the underlying y-array.
        """
        return self.y[0]

    def ymax(self):
        """Return the maximum of the underlying y-array.
        """
        return self.y[-1]

    def scale(self, scale_factor, zlabel=None):
        """Scale the spline z values.
        """
        _x = numpy.copy(self.x)
        _y = numpy.copy(self.y)
        _z = numpy.copy(self.z) * scale_factor
        if zlabel is None:
            zlabel = self.zlabel
        kx, ky = self.degrees
        return self.__class__(_x, _y, _z, None, kx, ky, self.xlabel, self.ylabel, zlabel)

    def __call__(self, x, y, dx=0, dy=0, grid=False):
        """Overloaded __call__method.

        Here we basically override the default value of the `grid` parameter
        from `True` to `False`, since we're typically interested in evaluating
        the spline at given physical coordinates, rather than grid points.
        """
        return RectBivariateSpline.__call__(self, x, y, dx=dx, dy=dy, grid=grid)

    def vslice(self, x, k=3):
        """Return a vertical slice at a given x of the bivariate spline.

        Args
        ----
        x : float
            The x value at which the vertical slice should be calculated.

        k : int, optional
            The degree of the resulting spline.
        """
        assert x >= self.xmin() and x <= self.xmax()
        _x = self.y
        _y = self(x, _x)
        fmt = dict(xlabel=self.ylabel, ylabel=self.zlabel)
        return xInterpolatedUnivariateSpline(_x, _y, k=k, **fmt)

    def hslice(self, y, k=3):
        """Return an horizontal slice at a given y of the bivariate spline.

        Args
        ----
        y : float
            The y value at which the horizontal slice should be calculated.

        k : int, optional
            The degree of the resulting spline.
        """
        assert y >= self.ymin() and y <= self.ymax()
        _x = self.x
        _y = self(_x, y)
        fmt = dict(xlabel=self.xlabel, ylabel=self.zlabel)
        return xInterpolatedUnivariateSpline(_x, _y, k=k, **fmt)

    def build_horizontal_ppf(self):
        """Create the horizontal percent point function (or inverse of cdf).

        Warning
        -------
        This really, really need to be fixed. Instead of grabbing a vertical
        slice at xmean, we should pass an argument to the function so that
        the subclasses can implement whatever is right for them.
        """
        _ymean = 0.5 * (self.ymin() + self.ymax())
        _refppf = self.hslice(_ymean).build_ppf()
        _x = _refppf.x
        _y = self.y.copy()
        _z = numpy.zeros(shape=(_x.size, _y.size))
        for j, _yp in enumerate(_y):
            _ppf = self.hslice(_yp).build_ppf()
            for i, _xp in enumerate(_x):
                _z[i, j] = _ppf(_xp)
        fmt = dict(xlabel='q', ylabel=self.ylabel, zlabel=self.zlabel)
        return xInterpolatedBivariateSplineLinear(_x, _y, _z, **fmt)

    def plot_contours(self, num_contours=10, colors='black', fontsize='small',
        logz=False, cfmt='%1.3f'):
        """Contour plot.
        """
        x, y = self.transposed_meshgrid(self.x, self.y)
        kwargs = dict(colors=colors)
        if logz:
            kwargs['norm'] = matplotlib.colors.LogNorm()
        cont = plt.contour(x, y, self.z, num_contours, **kwargs)
        plt.gca().clabel(cont, fontsize=fontsize, fmt=cfmt)

    def plot_color(self, num_contours=75, logz=False, **kwargs):
        """Color plot
        """
        x, y = self.transposed_meshgrid(self.x, self.y)
        vmin = kwargs.get('vmin', self.z.min())
        vmax = kwargs.get('vmax', self.z.max())
        if logz:
            logmin = numpy.log10(vmin)
            logmax = numpy.log10(vmax)
            levels = numpy.logspace(logmin, logmax, num_contours)
            cont_plot = plt.contourf(x, y, self.z, levels=levels,
                                     norm=matplotlib.colors.LogNorm())
        else:
            cont_plot = plt.contourf(x, y, self.z, num_contours, **kwargs)
        cont_plot.set_array(self.z)
        cont_plot.set_clim(vmin, vmax)
        if logz:
            cont_plot.set_norm(matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))
        color_bar = plt.colorbar(cont_plot, ax=plt.gca())
        color_bar.ax.zorder = -1
        if self.zlabel is not None:
            color_bar.set_label(self.zlabel)

    def plot(self, num_contours=75, logz=False, **kwargs):
        """Plot the spline.
        """
        self.plot_color(num_contours, logz, **kwargs)
        setup_gca(self.xlabel, self.ylabel)



class xInterpolatedBivariateSplineLinear(xInterpolatedBivariateSpline):

    """Bivariate linear interpolated spline on a rectangular grid.
    """

    def __init__(self, x, y, z, xlabel=None, ylabel=None, zlabel=None):
        fmt = dict(xlabel=xlabel, ylabel=ylabel, zlabel=zlabel)
        xInterpolatedBivariateSpline.__init__(self, x, y, z, None, 1, 1, **fmt)



class xStepFunction:

    """Small convenience class describing a step function.

    The basic rule, here, is that the value of the function between x[i] and x[i + 1]
    is exactly y[i].

    Args
    ----
    x : array_like
        The array of x values.

    y : array_like
        The array of y values, with shape (N, ), assuming that (N + 1, ) is the x shape.

    xlabel : str
        The text label for the x axis.

    ylabel : str
        The text label for the y axis.
    """

    def __init__(self, x, y, xlabel=None, ylabel=None):
        """Constructor.
        """
        assert len(x) == len(y) + 1
        self.x = x
        self.y = y
        self.xlabel = xlabel
        self.ylabel = ylabel

    def __call__(self, x):
        """Overloaded method.
        """
        if isinstance(x, Number):
            x = numpy.array(x)
        value = numpy.zeros(x.shape)
        index = (numpy.searchsorted(self.x, x) - 1).clip(0)
        for i in range(len(self.y)):
            value[index==i] = self.y[i]
        return value

    def plot(self, annotate=False, fmt='%.2f'):
        """Plot the step function.

        See the interesting discussion of a related problem at
        https://stackoverflow.com/questions/5347065
        (For the reference, this is Will's solution.)
        """
        x = numpy.ravel((self.x, self.x), order='F')[1:-1]
        y = numpy.ravel((self.y, self.y), order='F')
        plt.plot(x, y)
        setup_gca(xlabel=self.xlabel, ylabel=self.ylabel)
        ymin, ymax = plt.ylim()
        delta = (ymax - ymin) * 0.005
        if annotate:
            for _x, _y in zip(0.5 * (self.x[:-1] + self.x[1:]), self.y):
                plt.text(_x, _y + delta, fmt % _y, ha='center', color=last_line_color())

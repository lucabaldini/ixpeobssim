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

"""Custom random generation classes, mostly based on splines.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.core.spline import xInterpolatedBivariateSpline
from ixpeobssim.utils.logging_ import logger, abort

# pylint: disable=invalid-name, too-many-arguments

class xUnivariateGenerator(xInterpolatedUnivariateSpline):

    """Univariate random number generator based on a linear interpolated
    spline.

    This class is adding very little on top of the xInterpolatedUnivariateSpline
    spline class: essentially we are building the ppf in the constructor and
    we are adding the pdf() method (simply an alias to the underlying __call__()
    operator, as well as an rvs() method to generate actual random numbers).

    Args
    ----
    rv : (N,) array_like
        Array of points sampling the values of the random variable.

    pdf : (N,) array_like
        pdf values at the array rv.

    w : (N,) array_like, optional
        Weights for spline fitting.

    bbox : (2,) array_like, optional
        Bounduary of the approximation interval. .

    k : int
        The order of the spline (<=5, default 3, i.e., a cubic spline).

    ext : int or string, optional
        Controls the extrapolation mode for elements not in the interval defined
        by the knot sequence.

    xlabel: str, optional
        A text label of the random variable.

    ylabel: str, optional
        A text label for the pdf.
    """

    def __init__(self, rv, pdf, w=None, bbox=None, k=3, ext=0, xlabel='rv',
                 ylabel='pdf'):
        """ Constructor.
        """
        args = (rv, pdf, w, bbox, k, ext, xlabel, ylabel)
        xInterpolatedUnivariateSpline.__init__(self, *args)
        # This is a pdf, so we want to ensure that the underlying y values are
        # all positive (semi-definite).
        mask = self.y < 0.
        if mask.sum() > 0:
            logger.error('Passing %d negative values (out of %d) to %s',
                         mask.sum(), self.y.size, self.__class__.__name__)
            idx = numpy.argwhere(mask)
            logger.error('%s: %s', xlabel, self.x[idx].T)
            abort('A probability density function must be positive-defined')
        self.cdf = self.build_cdf()
        self.ppf = self.build_ppf()

    def rvs(self, size=1):
        """Return random variates of arbitrary size.

        Args
        ----
        size : int
            The number of random numbers to be generated.
        """
        return self.ppf(numpy.random.sample(size))

    def rvs_bounded(self, size=1, rvmin=None, rvmax=None):
        """Return random variates of arbitrary size between the specified bounds.

        This is an alternative to rvs(), where the uniformly-distributed random
        array passed to the ppf is spanning the appropriate subset of the [0, 1]
        interval.

        Although the default behavior of this method is identical to the standard
        rvs(), we prefer to have a separate function to avoid unnecessary
        complexity in what is by far the most common use case.
        """
        if rvmin is None:
            umin = 0.
        else:
            umin = self.cdf(rvmin)
        if rvmax is None:
            umax = 1.
        else:
            umax = self.cdf(rvmax)
        u = numpy.random.uniform(umin, umax, size)
        return self.ppf(u)



class xUnivariateGeneratorLinear(xUnivariateGenerator):

    """Subclass of xUnivariateGenerator restricted to a linear spline.
    """

    def __init__(self, rv, pdf, ext=0, xlabel='rv', ylabel='pdf'):
        """ Constructor.
        """
        args = (rv, pdf, None, None, 1, ext, xlabel, ylabel)
        xUnivariateGenerator.__init__(self, *args)



class xUnivariateAuxGenerator(xInterpolatedBivariateSpline):

    """Univariate random generator where the pdf of the random variable
    depends on an additional auxiliary variable.

    Internally, a proper meshgrid of random and auxiliary variable values is
    created, and the pdf is evaluated on the grid to construct an interpolated
    bivariate spline, so that a horizontal slice of the bivariate spline at any
    given value of the auxiliary variable is the actual pdf of the random
    variable at that value of the auxiliary variable.

    Note that the z function argument can either be an array-like object with
    shape (rv.size, aux.size), or a function that can be evaluated in the
    bi-dimensional phase space---see the base class xInterpolatedBivariateSpline
    for the implementation details.

    Args
    ----
    rv : array_like
        Input values for the random variable (assumed to be sorted).

    aux : array_like
        Input values for the auxiliary variable (assumed to be sorted).

    pdf : array_like of shape (x.size, y.size) or callable
        Input pdf values, with shape (x.size,y.size).

    bbox : array_like, optional
        The boundary of the approximation domain.

    kx, ky : int
        The spline order on the two axes.

    xlabel: str, optional
        A text label for the random variable.

    ylabel: str, optional
        A text label for the auxiliary variable.

    zlabel: str, optional
        A text label for the pdf values.

    Note
    ----
    A major backward-incompatible change was introduced here in version 3.0.0,
    in the context of https://bitbucket.org/ixpesw/ixpeobssim/issues/196,
    and the constructor signature of the class now reads as
    xUnivariateAuxGenerator(rv, aux, pdf...), instead of
    xUnivariateAuxGenerator(aux, rv, pdf...) as it used to be historically in
    all the previous versions of ixpeobssim.

    (On a related note, since some of the arguments are passed by name in
    the base spline classes---mostly the axis labels---we decided to keep the
    argument names in synch with those of the base classes.)

    The main driver for making this breaking change was that one the main use
    cases for this class is to generate time-dependent spectra, and the basic
    signature that we picked in this case, as specified in the documentation,
    is spec(E, t)---where E (the energy) is the random variable and t (the time)
    the auxiliary variable. This makes sense, as in most cases the spectra are
    time-independent, and keeping the energy in the first position allows to
    give the time a default value, simplifying all the signatures in this fairly
    common use case. This also makes the signature of the pdf(rv, aux) class
    method more reasonable.

    After the fix to the creation of the underlying meshgrid implemented in
    https://bitbucket.org/ixpesw/ixpeobssim/commits/cdbd99a7545bb1b0fd09591b2c9b973fb56325d4
    insisting on the signature as xUnivariateAuxGenerator(aux, rv, pdf...)
    would have forced us to swap the spectrum arguments spec(E, t) -> spec(t, E)
    all over the place (e.g., in srcmodel.spectrum).

    (In brief: we have rotated the underlying spline representation by 90
    degrees with respect to the historical implementation.)
    """
    def __init__(self, rv, aux, pdf, bbox=None, kx=3, ky=3,
                 xlabel='rv', ylabel='aux', zlabel='pdf'):
        """Constructor.
        """
        args = (rv, aux, pdf, bbox, kx, ky, xlabel, ylabel, zlabel)
        xInterpolatedBivariateSpline.__init__(self, *args)
        # This is a pdf, so we want to ensure that the underlying z values are
        # all positive (semi-definite).
        mask = self.z < 0.
        if mask.sum() > 0:
            logger.error('Passing %d negative values (out of %d) to %s',
                         mask.sum(), self.z.size, self.__class__.__name__)
            xidx , yidx = numpy.hsplit(numpy.argwhere(mask), 2)
            logger.error('%s: %s', xlabel, self.x[xidx].T)
            logger.error('%s: %s', ylabel, self.y[yidx].T)
            abort('A probability density function must be positive-defined')
        self.ppf = self.build_horizontal_ppf()

    def slice(self, aux):
        """Return the one-dimensional pdf for a given value of the auxiliary
        variable (i.e., a horizontal slice of the underlying bivariate spline).
        """
        return self.hslice(aux)

    def rvs(self, aux):
        """Return random variates for a given array of values of the auxiliary
        variable.

        Note that the sample size, here, is determined by the size of the
        aux argument.
        """
        return self.ppf(numpy.random.sample(len(aux)), aux)

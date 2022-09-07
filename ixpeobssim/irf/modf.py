# Copyright (C) 2015--2022, the ixpeobssim team.
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

"""Modulation-related facilities.
"""

from __future__ import print_function, division

import numbers

import numpy

from ixpeobssim.core.rand import xUnivariateAuxGenerator
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline, xInterpolatedBivariateSpline
from ixpeobssim.irf.base import xSpecRespBase
from ixpeobssim.utils.logging_ import logger, abort


# pylint: disable=invalid-name, too-many-ancestors, no-member


class xAzimuthalResponseGenerator(xUnivariateAuxGenerator):

    """Random number generator for the azimuthal response of the polarimeter.

    Here is the basic underlying math. Typically the response of a polarimeter
    to monochromatic, linearly polarized incident radiation is of the form:

    .. math:: N(\\phi) = A + B \\cos^2(\\phi - \\phi_0).

    This can be conveniently rewritten in term of the overall normalization
    (i.e., the total number of events) and the modulation, defined as

    .. math::
        m = \\frac{N_\\text{max} - N_\\text{min}}
        {N_\\text{max} + N_\\text{min}} = \\frac{B}{2A + B}

    (the modulation can also be characterized as the fraction of modulated
    events in a given distribution, as can be readily demonstrated, and it is
    the quantity that the detector is effectively measuring).
    The angular response can be rewritten as

    .. math::
        N(\\phi) = \\frac{N_0}{2\\pi} \\left[
        1 + m \\cos(2(\\phi - \\phi_0))
        \\right].

    For completeness, the modulation, the modulation factor and the polarization
    degree for a monocromatic source are related to each other through:

    .. math::
        m(E, t) = P(E, t) \\times \\mu(E)

    (i.e., at any given energy the modulation factor is the modulation of the
    detector response for 100% linearly polarized incident radiation).

    In terms of throwing random numbers, the phase is a trivial constant that
    can be added after the fact (modulo 2pi), so effectively the
    relevant probability density function is

    .. math::
        \\text{pdf}(\\phi) = \\frac{1}{2\\pi} \\left[
        1 + m \\cos(2\\phi) \\right],

    .. image:: ../docs/figures/test_azimuthal_resp_pdf.png

    The corresponding cumulative distribution function is

    .. math::
        \\text{cdf}(\\phi) = \\frac{1}{2\\pi} \\left(
        \\phi + \\frac{m}{2}\\sin{(2\\phi)} \\right),

    and it is unfortunate that it cannot be inverted, otherwise we would
    have no need to interpolate for generating random numbers according to
    this distribution.

    .. image:: ../docs/figures/test_azimuthal_resp_cdf.png

    From the prospecive of the code, this generator is a standard
    `xUnivariateAuxGenerator` where the azimuthal angle is our
    random variable and the modulation is our auxiliary variable. For any given
    value of the modulation, a vertical slice is providing the corresponding
    one-dimensional pdf.

    .. image:: ../docs/figures/test_azimuthal_resp_generator.png

    The class also provide facilities to fit a histogram to recover the
    underlying modulation and phase.

    Example
    -------
    >>> import numpy
    >>> from ixpeobssim.utils.matplotlib_ import plt
    >>> from ixpeobssim.irf.modf import xAzimuthalResponseGenerator
    >>>
    >>> generator = xAzimuthalResponseGenerator()
    >>> modulation = numpy.full(1000000, 0.5)
    >>> phase = numpy.radians(45.)
    >>> phi = generator.rvs_phi(phase, modulation)
    >>> hist = plt.hist(phi, bins=numpy.linspace(0, 2*numpy.pi, 100),
                        histtype='step', label='Random angles')
    >>> fit_results = generator.fit_histogram(hist)
    >>> fit_results.plot()
    >>> plt.show()

    .. image:: ../docs/figures/test_azimuthal_resp_rvs.png
    """

    def __init__(self):
        """Constructor.
        """
        phi = numpy.linspace(-numpy.pi, numpy.pi, 100)
        m = numpy.linspace(0., 1., 100)
        fmt = dict(xlabel='$\\phi$ [rad]', ylabel='$m$')
        xUnivariateAuxGenerator.__init__(self, phi, m, self.pdf, **fmt)

    @staticmethod
    def pdf(phi, m):
        """Evaluate the underlying one-dimensional pdf for a given value of the
        modulation, and assuming that the phase is identically zero.

        Arguments
        ---------
        phi : float or array
            The (independent) azimuthal angle variable, in [-pi, pi].

        m : float or array
            The modulation, in [0, 1].
        """
        return (1 + m * numpy.cos(2. * phi)) / (2. * numpy.pi)

    @staticmethod
    def cdf(phi, m):
        """Evaluate the underlying one-dimensional cdf for a given value of the
        modulation, and assuming that the phase is zero.

        Arguments
        ---------
        phi : float or array
            The (independent) azimuthal angle variable, in [-pi, pi].

        m : float or array
            The modulation, in [0, 1].
        """
        return 0.5 + (phi + 0.5 * m * numpy.sin(2. * phi)) / (2. * numpy.pi)

    def build_horizontal_ppf(self):
        """Overloaded method.

        This is essentially done for a few reasons:

        * we do have an analytical form of the cdf that we can use to
          reduce the numerical noise;
        * this random generator is peculiar in that the pdf is linear in the
          auxiliary variable m;
        * this is possibly the most important random engine in the framework
          and its accuracy must be tested carefully.

        There are many parameters that can be optimized, here, in order to
        maximize the accuracy of the ppf representation. (This accuracy is
        characterized in test/test_azimuthal_response.py in terms of standard
        and maximum relative error). Among these are the grids sampling the
        modulation and quantile values, and the order of the various splines
        involved. We put some thoughts into getting this right, but we cannot
        exclude that this can be improved.

        For the modulation values we just take a constant-pitch grid, while for
        the quantiles we calculate the cdf (for m = 1, which is where the
        deviations from a straight lines are largest) on a constant-pitch grid
        in angle.

        We empirically found that k=3 is the best spline index for the ppf
        at any given modulation value, while we get the maximum accuracy for
        kx = ky = 1 in the actual bivariate spline representing the ppf.
        """
        phi = numpy.linspace(-numpy.pi, numpy.pi, 200)
        m = numpy.linspace(0., 1., 200)
        q = self.cdf(phi, m=1.)
        z = numpy.zeros(shape = (q.size, m.size))
        for i, _m in enumerate(m):
            _ppf = xInterpolatedUnivariateSpline(self.cdf(phi, _m), phi, k=3)
            z[:, i] = _ppf(q)
        fmt = dict(xlabel=self.ylabel, ylabel='q', zlabel=self.xlabel)
        return xInterpolatedBivariateSpline(q, m, z, kx=1, ky=1, **fmt)

    def rvs_phi(self, m, phase):
        """Generate random variates for any modulation and phase values.

        This is essentially calling the underlying xUnivariateAuxGenerator.rvs()
        method passing the modulation array as an argument and adding the phase
        array (modulo 2pi) to the output.

        Arguments
        ---------
        m : array
            An array of modulation values. (The function returns an equal-length
            array of phi values.)

        phase : float or array
            The phase of the modulation. (This can either be a vector or an
            array of the same length as `modulation`.)
        """
        return numpy.mod(self.rvs(m) + phase, 2 * numpy.pi) - numpy.pi



class xModulationFactor(xSpecRespBase):

    """Class describing the modulation factor.

    The effective area is essentially a linear spline, with built-in facilities
    for evaluation and plotting.

    To zero-th order, an `xModulationFactor` instance is an object capable of
    evaluating itself in a point or in a grid of points, and capable of
    plotting itself.

    More interestingly, it can generate random `phi` values, given a vector
    of event energies and corresponding vectors (or simple floats) representing
    the polarization degree and angle corresponding to the energies themselves.
    Internally, any `xModulationFactor` has an `xAzimuthalResponseGenerator`
    object and when the `xModulationFactor.rvs_phi()` method is called,
    the polarization degree is multiplied by the modulation factor of the
    detector, evaluated at the right energy, and converted into a modulation
    value, after which the underlying `xAzimuthalResponseGenerator.rvs_phi()`
    is called.
    """

    Y_UNITS = ''
    Y_LABEL = 'Modulation factor'

    def __init__(self, file_path, extension='fits', k=1):
        """Overloaded constructor
        """
        xSpecRespBase.__init__(self, file_path, extension, k)
        self.generator = xAzimuthalResponseGenerator()

    def rvs_phi(self, energy, polarization_degree, polarization_angle):
        """Return random variates for a given array of values of energy,
        polarization degree and polarization angle.

        Arguments
        ---------
        energy : array
            An array of energy values. (The function returns an equal-length
            array of phi values.)

        polarization_degree : array or float
            The polarization degree, in [0--1]. (This can either be a vector
            or an array of the same length as `energy`.)

        polarization_angle : array or float
            The polarization angle, in radians. (This can either be a vector or
            an array of the same length as `energy`.)
        """
        # If the input polarization degree is a number (i.e., a float) convert
        # it to a numpy array of the same shape of the energy with identical
        # values. This makes our life easier downstream in verifying that the
        # input values for the polarization degrees lie between 0 and 1---and
        # to provide a sensible diagnostic message if this is not the case.
        if isinstance(polarization_degree, numbers.Number):
            polarization_degree = numpy.full(energy.shape, polarization_degree)
        # Make sure that the polarization degree is within the physical range.
        # Note that this chack is comparatively hard to do at the model level,
        # before the actual event energies have been extracted, as in principle
        # one would need to evaluate the model itself on an infinitely fine
        # grid.
        def _msg(condition):
            """Small nested function to print out some useful diagnostics about
            possible unphysical values for the polarization degree.
            """
            mask = numpy.where(condition)
            efunky = energy[mask]
            pfunky = polarization_degree[mask]
            msg = 'Unphysical polarization degree %s at %d input energies %s'
            logger.error(msg, pfunky, len(efunky), efunky)
        # Do the actual checks.
        if polarization_degree.min() < 0.:
            _msg(polarization_degree < 0.)
            abort('The polarization degree must be >= 0')
        if polarization_degree.max() > 1.:
            _msg(polarization_degree > 1.)
            abort('The polarization degree must be <= 1')
        # And now we're good to go.
        modulation = self(energy) * polarization_degree
        return self.generator.rvs_phi(modulation, polarization_angle)

    def weighted_average(self, energy):
        """Return the weighted average of the mudulation factor given an
        array of energies.

        Arguments
        ---------
        energy : array
            An array of energy values.
        """
        return self(energy).sum() / len(energy)

    def plot(self):
        """Overloaded method for the xpirfview application.
        """
        # pylint: disable=arguments-differ
        self.plot_base()

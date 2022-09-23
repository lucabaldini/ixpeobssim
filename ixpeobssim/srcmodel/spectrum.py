#!/urs/bin/env python
#
# Copyright (C) 2015--2020, the ixpeobssim team.
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

"""Spectral facilities.
"""

from __future__ import print_function, division

import sys

import numpy

from ixpeobssim.utils.environment import PYXSPEC_INSTALLED
from ixpeobssim.core.rand import xUnivariateGenerator, xUnivariateAuxGenerator
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.evt.mdp import xMDPTable
from ixpeobssim.utils.units_ import erg_to_keV, keV_to_erg
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.srcmodel.gabs import xInterstellarAbsorptionModel
from ixpeobssim.utils.fmtaxis import fmtaxis
import ixpeobssim.irf.ebounds as ebounds
from ixpeobssim.srcmodel import load_tabular_data
from ixpeobssim.utils.misc import pairwise
if PYXSPEC_INSTALLED:
    from ixpeobssim.evt.xspec_ import sample_spectral_model

# pylint: disable=invalid-name, too-many-arguments


def load_spectral_spline(file_path, emin=ebounds.ENERGY_MIN, emax=ebounds.ENERGY_MAX,
                         energy_col=0, flux_col=1, delimiter=None, **kwargs):
    """Convenience function to load spectral data from a txt file.

    Since we happen to load spectral data from text files a lot, this is an
    attempt to provide a facility that is generic enough to be effectively
    reused.
    """
    data = load_tabular_data(file_path, emin, emax, energy_col, delimiter)
    energy = data[energy_col]
    flux = data[flux_col]
    kwargs.update(fmtaxis.spec)
    return xInterpolatedUnivariateSpline(energy, flux, **kwargs)


def wrap_spectral_parameter(parameter):
    """Wrap a spectral parameter in order to handle time- (or phase-) dependence.

    This is a small, handy trick to handle in a consistent fashion spectral
    models where the underlying parameters can be either constant or
    time-dependent: once they are wrapped, the parameters can be called with a
    given time as an argument even when they are time-independent.

    Note that, in all cases, the t argument defaults to None, to allow
    time-independent spectral models to be called omitting it.

    Arguments
    ---------
    parameter : float or callable
        The spectral parameter (e.g., the normalization or the index of a
        power law). This can be either a scalar or a callable with the
        signature parameter(t).

    Returns
    -------
        An anonymous function that, called with a given time as the only
        arguments, returns the parameter evaluated at that time (or the scalar
        if the parameter is time-independent).

    Example
    -------
    >>> from ixpeobssim.srcmodel.spectrum import wrap_spectral_parameter
    >>>
    >>> param = wrap_spectral_parameter(1.)
    >>> print(param(1000.))
    >>> 1.0
    >>> print(param())
    >>> 1.0
    """
    if callable(parameter):
        return lambda t=None: parameter(t)
    return lambda t=None: parameter


def wrap_spectral_model(model):
    """Simple decorator to simplify the definition of spectral models.

    This is essentially wrapping all the model parameters with the
    wrap_spectral_parameter() function and calling the input models with the
    wrapped parameters, so that we don't have to worry about whether each of
    them is time-dependent or time-independent.
    """
    def wrapper(*args):
        """Wrapper definition.

        Warning
        -------
        We might use some functools magic, here, to make sure the decorated
        model factories preserve the correct function signature when queried for
        help.
        """
        args = (wrap_spectral_parameter(arg) for arg in args)
        return model(*args)
    return wrapper


@wrap_spectral_model
def power_law(norm, index):
    """Simple power law.
    """
    return lambda E, t=None: \
        norm(t) * numpy.power(E, -index(t))


@wrap_spectral_model
def cutoff_power_law(norm, index, cutoff):
    """Power law with a high-energy exponential cutoff.
    """
    return lambda E, t=None: \
        norm(t) * numpy.power(E, -index(t)) * numpy.exp(-E / cutoff(t))


@wrap_spectral_model
def highecut(ecut, efold):
    """highecut multiplicative component, see
    https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node244.html
    """
    return lambda E, t=None: \
        (E <= ecut(t)) + (E > ecut(t)) * numpy.exp((ecut(t) - E) / efold(t))


@wrap_spectral_model
def highecut_power_law(norm, index, ecut, efold):
    """Power law modulated with a highecut multiplicative component, see
    https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node244.html
    """
    return lambda E, t=None: \
        power_law(norm(t), index(t))(E, t) * highecut(ecut(t), efold(t))(E, t)


def gaussian_line(norm, mean, sigma):
    """Gaussian line spectral component.
    """
    return lambda E, t=None: \
        norm / numpy.sqrt(2. * numpy.pi) / sigma * numpy.exp(-(E - mean)**2. / 2. / sigma**2.)


# Aliases for compatibility with XSPEC, in case anybody cares.
powerlaw = power_law
cutoffpl = cutoff_power_law
gauss = gaussian_line




class xXspecModel(xInterpolatedUnivariateSpline):

    """Basic interface to a generic time-independent XSPEC model.

    This is essentially a univariate interpolated spline that is constructed
    from a regular-grid sample of a generic XSPEC model, defined by
    an expression and a (full) set of parameters.

    The model can be loaded from a text file in a suitable form using
    the xXspecModel.from_file() facility, and the data points can be exported
    to a text file using the save_ascii() class method.

    Arguments
    ---------
    expression : str
        The model expression string, using full XSPEC component names.

    parameters : dict or list
        The model parameters. This can take the form of either a list, where
        all the parameters are passed (in the right order), or a dictionary
        indexed by the serial identifier of the parameter itself.
        The second form allows for passing along only a subset of the parameters
        and mimics the behavior of the XSPEC Python bindings described in
        https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/model.html

        Note that any attempt of implementing a structure where the parameters
        are passed by name is made non trivial by the possibility of compound
        models featuring multiple instances of the same parameter names (for
        different components), and the questionable benefit provided by this
        idea was deemed not worth the additional complications connected with
        that.

    w, bbox, k, ext
        Set of parameters controlling the spline

    emin, emax : double
        Energy limits, in keV.

    num_points : int
        The number of points to sample the spectrum.

    correct : bool
        If True, we attempt at converting the integral values from XSPEC into
        actual fluxes at the bin centers. (This is achieved via a numerical
        integration of a cubic spline.)
    """

    def __init__(self, expression, parameters, w=None, bbox=None, k=3, ext=0,
                 emin=ebounds.ENERGY_MIN, emax=ebounds.ENERGY_MAX, num_points=250,
                 correct=True):
        """Constructor.
        """
        # Mind that getting rid of the spaces here reduces the chances of having
        # strings with spaces passed around as command-line arguments at later
        # stages.
        self.expression = expression.replace(' ', '')
        self.parameters = parameters
        args = expression, parameters, emin, emax, num_points
        binning, energy, flux, self.parameter_names = sample_spectral_model(*args)
        if correct:
            # Build a first spline with the data from XSPEC.
            spline = xInterpolatedUnivariateSpline(energy, flux, k=3, ext=0)
            # Overwrite the spline including the extrapolation of the previous
            # one over the entire [emin, emax] range---this is needed because
            # we shall perform a numerical integration, and scipy assumes that
            # the spline is identically zero outside its domain.
            spline = xInterpolatedUnivariateSpline(binning, spline(binning), k=3, ext=0)
            bin_width = numpy.diff(binning)
            # Calculate the conversion factor between the average flux in the bin
            # and the value at the bin center.
            scale = numpy.array([spline.integral(binning[i], binning[i + 1]) \
                for i in range(len(energy))]) / bin_width / spline(energy)
            flux /= scale

        args = energy, flux, w, bbox, k, ext
        xInterpolatedUnivariateSpline.__init__(self, *args, **fmtaxis.spec)

    def __call__(self, E, t=None):
        """Overload __call__ dunder.

        This is done so that the spline can be fed as a spectrum object into
        xpobssim directly, without any need of wrap the arguments.
        """
        return xInterpolatedUnivariateSpline.__call__(self, E)

    @classmethod
    def from_file(cls, file_path, w=None, bbox=None, k=3, ext=0,
                  emin=ebounds.ENERGY_MIN, emax=ebounds.ENERGY_MAX, num_points=250):
        """Create a class instance from file.

        The basic file format is the following. (Note that it is up to the
        user to make sure that the components are defined in the same order
        in which they first appear in the model expression string, and that
        the order of the parameters within each component is exactly the
        one that XSPEC is expecting.)

        ::

            model     phabs*powerlaw
            COMP1	  phabs
            nH	  3.
            COMP2     powerlaw
            PhoIndex  1.
            norm      100.
        """
        logger.info('Reading XSPEC model from %s...', file_path)
        expression = None
        parameters = []
        with open(file_path) as input_file:
            for line in input_file:
                line = line.strip()
                if line.startswith('#') or line.isspace() or len(line) == 0:
                    continue
                key, value = line.split(None, 1)
                if key.startswith('COMP'):
                    pass
                elif key == 'model':
                    expression = value
                else:
                    parameters.append(float(value))
        logger.info('Done, expression = "%s", parameters = %s', expression, parameters)
        assert expression is not None
        return cls(expression, parameters, w, bbox, k, ext, emin, emax, num_points)

    def save_ascii(self, file_path):
        """Save a txt file with spectrum data points.
        """
        logger.info('Saving XSPEC model to %s...', file_path)
        with open(file_path, 'w') as output_file:
            for e, f in zip(self.x, self.y):
                output_file.write('%.4e    %.4e\n' % (e, f))
        logger.info('Done.')

    def __str__(self):
        """String formatting.
        """
        text = 'Model "%s"' % self.expression
        for name, val in zip(self.parameter_names, self.parameters):
            text += '\n%s = %.3e' % (name, val)
        return text



def integral_flux(spectrum, emin, emax, column_density=None):
    """Return the integral flux within a given energy interval for a given
    spectral model.

    Arguments
    ---------
    spectrum : callable
        The spectral model (must be able to evaluate itself onto an array of energies)

    emin : float
        The minimum energy for the integral

    emax : float
        The minimum energy for the integral

    column_density : float (optional)
        The value of the column density (in cm^-2) used to calculate the
        Galactic absorption.
    """
    x = numpy.linspace(emin, emax, 1000)
    y = spectrum(x)
    if column_density is not None:
        y *= xInterstellarAbsorptionModel().transmission_factor(column_density)(x)
    return xInterpolatedUnivariateSpline(x, y).integral(emin, emax)


def integral_energy_flux(spectrum, emin, emax, column_density=None, erg=True):
    """Return the integral energy flux within a given energy interval for a
    given spectral model.

    Arguments
    ---------
    spectrum : callable
        The spectral model (must be able to evaluate itself onto an array of energies)

    emin : float
        The minimum energy for the integral

    emax : float
        The minimum energy for the integral

    column_density : float (optional)
        The value of the column density (in cm^-2) used to calculate the
        Galactic absorption.

    erg : bool
        If True, convert the output from keV to erg
    """
    value = integral_flux(lambda energy: energy * spectrum(energy), emin, emax, column_density)
    if erg:
        value = keV_to_erg(value)
    return value


def pl_integral(norm, index, emin, emax):
    """Return the integral of a generic power law in a given energy range.
    """
    return norm / (1. - index) * (emax**(1 - index) - emin**(1 - index))


def pl_integral_flux(norm, index, emin, emax):
    """Return the integral flux for a generic power law in a given energy range.
    """
    return keV_to_erg(pl_integral(norm, index - 1., emin, emax))


def pl_norm(integral, emin, emax, index, energy_power=0.):
    """Return the power-law normalization resulting in a given integral
    flux (or integral energy flux, or more in general integral of the
    flux multiplied by a generic power of the energy) between the minimum and
    maximum energies assuming a given spectral index.

    More specifically, given a power law of the form

    .. math::
       \\mathcal{S}(E) = C\\left( \\frac{E}{E_0} \\right)^{-\\Gamma}
       \\quad [{\\rm keV}^{-1}~{\\rm cm}^{-2}~{\\rm s}^{-1}],

    (where :math:`E_0 = 1~{\\rm keV}`) we define
    :math:`\\beta = (1 + p - \\Gamma)` and calculate

    .. math::
       I_{p} = \\int_{E_{\\rm min}}^{E_{\\rm max}} E^{p}\\mathcal{S}(E) dE =
       \\begin{cases}
       \\frac{C E_0^{\\Gamma}}{\\beta}
       \\left( E_{\\rm max}^{\\beta} - E_{\\rm min}^{\\beta}\\right)
       \\quad \\beta \\neq 0\\\\
       C E_0^{\\Gamma} \\ln \\left( E_{\\rm max}/E_{\\rm min} \\right)
       \\quad \\beta = 0\\\\
       \\end{cases}
       \\quad [{\\rm keV}^{p}~{\\rm cm}^{-2}~{\\rm s}^{-1}].

    Hence

    .. math::
        C =
        \\begin{cases}
        \\frac{I_p\\beta}{E_0^{\\Gamma}
        \\left( E_{\\rm max}^{\\beta} - E_{\\rm min}^{\\beta}\\right)}
        \\quad \\beta \\neq 0\\\\
        \\frac{I_p}{E_0^{\\Gamma}
        \\ln \\left( E_{\\rm max}/E_{\\rm min} \\right)}
        \\quad \\beta = 0.
        \\end{cases}

    Arguments
    ---------
    integral : float or array
        The value of the integral flux or energy-to-some-power flux

    emin : float
        The minimum energy for the integral flux

    emax : float
        The maximum energy for the integral flux

    index : float
        The power-law index

    energy_power : float
        The power of the energy in the integral

    erg : bool
        if True, convert erg to keV in the calculation.
    """
    assert emax > emin
    beta = 1 + energy_power - index
    if beta != 0:
        return integral * beta / (emax**beta - emin**beta)
    return integral / numpy.log(emax / emin)


def int_eflux2pl_norm(integral, emin, emax, index, erg=True):
    """Convert an integral energy flux into the corresponding power-law
    normalization.
    """
    if erg:
        integral = erg_to_keV(integral)
    return pl_norm(integral, emin, emax, index, 1.)



class xSourceSpectrum(xUnivariateAuxGenerator):

    """Class representing a source spectrum,

    .. math::
       \\mathcal{S}(E, t)
       \\quad [\\text{cm}^{-2}~\\text{s}^{-1}~\\text{keV}^{-1}].

    At the top level this is essentially a bivariate spline with the energy
    running on the x-axis, the time running on the y-axis, and the differential
    source spectrum on the x-axis.

    Arguments
    ---------
    E : array_like
        The energy array (in keV) where the spectrum is evaluated.

    t : array_like
        The time array (in s) where the spectrum is evaluated.

    spectrum : callable
        The source spectrum, i.e., a Python function with the signature
        spectrum(E, t) returning the differential energy spectrum of the source
        at the appropriate enenrgy and time value(s).

    column_density : float, defaults to 0.
        The H column density to calculate the insterstellar absorption.

    redshift : float, defaults to 0.
        The source redshift.

    kx : int (optional)
        The degree of the bivariate spline on the x axis.

    ky : int (optional)
        The degree of the bivariate spline on the y axis.
    """

    XLABEL = 'Energy [keV]'
    YLABEL = 'Time [s]'
    ZLABEL = 'dN/dE [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]'

    def __init__(self, E, t, spectrum, column_density=0., redshift=0., kx=3, ky=3):
        """Constructor.
        """
        args = E, t, self._pdf(spectrum, column_density, redshift), None, kx, ky
        fmt = dict(xlabel=self.XLABEL, ylabel=self.YLABEL, zlabel=self.ZLABEL)
        xUnivariateAuxGenerator.__init__(self, *args, **fmt)

    @staticmethod
    def _pdf(spectrum, column_density, redshift):
        """Private method returning the actual pdf to be used to create the
        array for the z-axis of the underlying bivariate spline.

        In this case the function is performing the convolution of the source
        spectrum with the column density and the application of the redshift, as
        appropriate.

        Note that this is taking a spectrum (E, t) callable as its first
        argument and returning a callable (more precisely, a lambda function)
        with the same signature in output.

        Subclasses can overload the method and modify the behavior (i.e., add
        the convolution with the effective area in order to generate a count
        spectrum).
        """
        pdf = lambda E, t: spectrum(E * (1. + redshift), t)
        if column_density <= 0.:
            return pdf
        ism = xInterstellarAbsorptionModel()
        trans = ism.transmission_factor(column_density)
        return lambda E, t: trans(E) * pdf(E, t)

    def emin(self):
        """Return the minimum energy for which the spectrum is defined.
        """
        return self.xmin()

    def emax(self):
        """Return the maximum energy for which the spectrum is defined.
        """
        return self.xmax()

    def tmin(self):
        """Return the minimum time for which the spectrum is defined.
        """
        return self.ymin()

    def tmax(self):
        """Return the maximum time for which the spectrum is defined.
        """
        return self.ymax()

    def time_slice(self, t):
        """Return a one-dimensional slice of the spectrum at a given time.
        """
        return self.hslice(t)

    def energy_slice(self, E):
        """Return a one-dimensional slice of the spectrum at a given energy.
        """
        return self.vslice(E)

    def build_time_integral(self, tmin=None, tmax=None):
        """Build the time-integrated spectrum between the give bounds.

        The output is stored in the form of a xUnivariateGenerator, featuring
        all the spline facilities, along with the capability of extracting
        random numbers.
        """
        if tmin is None:
            tmin = self.tmin()
        if tmax is None:
            tmax = self.tmax()
        x = self.x
        y = numpy.array([self.energy_slice(E).integral(tmin, tmax) for E in x])
        ylabel = 'Time-integrated spectrum (%.1f--%.1f s)' % (tmin, tmax)
        fmt = dict(xlabel=self.xlabel, ylabel=ylabel)
        return xUnivariateGenerator(x, y, **fmt)

    def build_time_average(self, tmin=None, tmax=None):
        """Build the time-averaged spectrum.
        """
        if tmin is None:
            tmin = self.tmin()
        if tmax is None:
            tmax = self.tmax()
        x = self.x
        y = numpy.array([self.energy_slice(E).integral(tmin, tmax) for E in x])
        y /= (tmax - tmin)
        ylabel = 'Time-averaged spectrum (%.1f--%.1f s)' % (tmin, tmax)
        fmt = dict(xlabel=self.xlabel, ylabel=ylabel)
        return xUnivariateGenerator(x, y, **fmt)

    def build_energy_integral(self, emin=None, emax=None):
        """Build the energy-integrated spectrum between the given bounds.

        The output is stored in the form of a xUnivariateGenerator, featuring
        all the spline facilities, along with the capability of extracting
        random numbers.
        """
        if emin is None:
            emin = self.emin()
        if emax is None:
            emax = self.emax()
        x = self.y
        y = numpy.array([self.time_slice(t).integral(emin, emax) for t in x])
        ylabel = 'Energy-integrated spectrum (%.2f--%.2f keV)' % (emin, emax)
        fmt = dict(xlabel=self.ylabel, ylabel=ylabel)
        return xUnivariateGenerator(x, y, **fmt)

    def build_light_curve(self, emin=None, emax=None):
        """Build the light curve, i.e., the spectrum integrated over the
        entire energy range, as a function of time.
        """
        return self.build_energy_integral(emin, emax)



class xCountSpectrum(xSourceSpectrum):

    """Class representing a count spectrum, i.e., the convolution of the
    source photon spectrum and the detector effective area

    .. math::
       \\mathcal{C}(E, t) = \\mathcal{S}(E, t) \\times A_{\\rm eff}(E)
       \\quad [\\text{s}^{-1}~\\text{keV}^{-1}].

    Arguments
    ---------
    spectrum : callable
        The source spectrum, i.e., a Python function that can be called
        with two numpy arrays E and t and returns the differential energy
        spectrum of the source at the appropriate enenrgy and time value(s).

    aeff : ixpeobssim.irf.xEffectiveArea instance
        The instrument effective area.

    t : array_like
        The array of time values where we're sampling the light curve of the
        source.

    column_density : float, defaults to 0.
        The H column density to calculate the insterstellar absorption.

    redshift : float, defaults to 0.
        The source redshift.

    scale_factor : float, defaults to 1.
        An optional scale factor for the output count spectrum (see the
        warning below).

    emin : float (optional)
        The minimum value for the energy grid.

    emax : float (optional)
        The maximum value for the energy grid.

    kx : int (optional)
        The degree of the bivariate spline on the x axis.

    ky : int (optional)
        The degree of the bivariate spline on the y axis.

    Warning
    -------
    The scale_factor parameter was a workaround for running the MDP script on
    a pulsar, where we sample in phase and then have to multiply by the
    observation time. The functionality should be implemented in a more roboust
    fashion.
    """

    ZLABEL = 'dN/dE [s$^{-1}$ keV$^{-1}$]'

    def __init__(self, spectrum, aeff, t, column_density=0., redshift=0.,
                 scale_factor=1., emin=None, emax=None, kx=3, ky=3):
        """Constructor.
        """
        logger.info('Creating the count spectrum...')
        self.scale_factor = scale_factor
        # Build the convolution of the source spectrum with the effective
        # area (and the global scale factor, until we get rid of it).
        conv = lambda E, t: scale_factor * aeff(E) * spectrum(E, t)
        # Setup the array for the energy axis. Originally this was simply
        # taking the x array of the effective area table, then we added two
        # optional command-line switches to allow for running simulation in a
        # custom (restricted) energy range---see
        # https://bitbucket.org/ixpesw/ixpeobssim/issues/267
        if emin is None:
            emin = aeff.xmin()
        if emax is None:
            emax = aeff.xmax()
        # Note that at this point the number of points for the energy grid is
        # completly arbitrary. 1000 is quite generous, but when we moved away
        # from the original effective area array we got a couple of unit test
        # failures in correspondence of the Cu edge, and since the complexity
        # of a spline evaluation is logarithmic in the number of points we
        # decided to play safe.
        E = numpy.linspace(emin, emax, 1000)
        # Initialize the base class. Note that we intercept a possible
        # SystemExit exception from the parent class initialization in order
        # to provide a useful error message to the user, should something
        # go wrong here.
        try:
            xSourceSpectrum.__init__(self, E, t, conv, column_density, redshift, kx, ky)
        except SystemExit as e:
            print(e)
            msg = 'It looks like your input spectral model evaluates to '\
                  'negative values somewhere\nin the %.2f--%.2f keV energy band '\
                  '(the previous message might offer a clue as\nto where this is '\
                  'happening).\nPlease ensure that your model is '\
                  'positive-definite, or override the energy\nbounds for the '\
                  'simulation via the --emin and/or --emax command-line '\
                  'switches.' % (emin, emax)
            sys.exit(msg)
        # Build the light curve (we *always* need it).
        self.light_curve = self.build_light_curve()

    def rvs_event_times(self, size):
        """Extract a series of random event times from the count spectrum.

        Note the array is sorted in place.
        """
        time_ = self.light_curve.rvs(size)
        time_.sort()
        return time_

    def num_expected_counts(self, tmin=None, tmax=None, emin=None, emax=None):
        """Return the number of expected counts within a given time interval
        and energy range

        .. math::
           \\int_{t_{\\rm min}}^{t_{\\rm max}}
           \\int_{E_{\\rm min}}^{E_{\\rm max}}
           \\mathcal{C}(E, t) dt dE
        """
        if emin is None:
            emin = self.emin()
        if emax is None:
            emax = self.emax()
        if tmin is None:
            tmin = self.tmin()
        if tmax is None:
            tmax = self.tmax()
        return self.integral(emin, emax, tmin, tmax)

    def build_mdp_table(self, ebinning, modulation_factor):
        """Calculate the MDP values in energy bins, given the modulation
        factor of the instrument as a function of the energy.

        Arguments
        ---------
        ebinning : array
            The energy binning

        modulation_factor : ixpeobssim.irf.modf.xModulationFactor instance
            The instrument modulation factor as a function of the energy.
        """
        # Build the time-integrated spectrum.
        time_integrated_spectrum = self.build_time_integral()
        # Build the modulation-factor spectrum, i.e. the object that we
        # integrate to calculate the effective modulation factor in a
        # given energy bin for a given spectrum.
        _x = time_integrated_spectrum.x
        _y = time_integrated_spectrum.y * modulation_factor(_x)
        mu_spectrum = xInterpolatedUnivariateSpline(_x, _y, k=1)
        # Loop over the energy bins and calculate the actual MDP values.
        # Note that we also calculate the MDP on the entire energy range.
        observation_time = self.scale_factor * (self.tmax() - self.tmin())
        mdp_table = xMDPTable(observation_time)
        for emin, emax in pairwise(ebinning):
            num_counts = self.num_expected_counts(emin=emin, emax=emax)
            mu_eff = mu_spectrum.integral(emin, emax) / num_counts
            mdp_table.add_row(emin, emax, mu_eff, num_counts)
        return mdp_table

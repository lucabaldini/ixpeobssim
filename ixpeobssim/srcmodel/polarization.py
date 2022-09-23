#!/urs/bin/env python
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

"""Facilities related to modeling of the polarization properties.
"""

from __future__ import print_function, division

import numbers

import numpy
from astropy.io import fits
from astropy import wcs
from scipy.interpolate import RegularGridInterpolator

from ixpeobssim.utils.environment import ASTROPY_VERSION
from ixpeobssim.core.spline import xInterpolatedBivariateSplineLinear
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.srcmodel.img import xFITSImage
from ixpeobssim.utils.matplotlib_ import plt, plot_arrows
from ixpeobssim.core.stokes import xModelStokesParameters
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.units_ import arcmin_to_degrees, degrees_to_arcmin
from ixpeobssim.utils.math_ import modulo_2pi
from ixpeobssim.utils.astro import square_sky_grid



# pylint: disable=invalid-name, too-many-locals, too-many-arguments

INTERPOLATION_POS_BOUND_ERROR_MSG = \
"""%d of the input sky coordinates (%.3f%%) happen to be outside the Stokes
sky-map underlying interpolation grid (this might indicate that the FITS images
used for the polarization model are too small). Input positions outside the
bounds will be clipped to the edge of the interpolation grid, but you should
definitely consider looking into this."""

INTERPOLATION_ENE_BOUND_ERROR_MSG = \
"""%d of the input event energies (%.3f%%) happen to be oustide the Stokes
sky-cube energy layers. Offending energies will be clipped to the cube bounds,
but you should definitely consider extending the latters."""


def constant(C):
    """Simple wrapper returning a constant, independently of the input
    arguments.

    Args
    ----
    C : int or float
        The target constant.

    New in version x.x.x: we added default values to all the arguments in order
    to be able to call the function with no argument at all. In addition, when
    the first argument (i.e., the energy) is a numpy array we are now returning
    an array of the same length (whose elements are all equal to the target
    constant).
    """
    def _function(E=None, t=None, ra=None, dec=None):
        if isinstance(E, numpy.ndarray):
            return numpy.full(E.shape, C)
        if isinstance(t, numpy.ndarray):
            return numpy.full(t.shape, C)
        return C
    return _function


def constant_spline(xmin, xmax, C, xlabel='Energy [keV]', ylabel=None):
    """Convenience function to generate a constant spline on a generic
    x-axis interval. (Note that underlying array only have two points, as
    in this case the interpolation is trivial.)

    Args
    ----
    x : numpy array
        The grid for the x-axis of the output spline.

    C : int or float
        The target constant

    We typically use this to express a polarization degree or angle as a
    function of the energy, which is the reason for the default argument of
    the name and units on the x-axis.
    """
    x = numpy.array([xmin, xmax])
    y = numpy.full(x.shape, C)
    return xInterpolatedUnivariateSplineLinear(x, y, xlabel=xlabel, ylabel=ylabel)



def _broadband_calc_base(spec, pol_func, emin=2., emax=8., num_points=500):
    """Basic function for the calculation of broadband averages of polarization
    parameters.
    """
    E = numpy.linspace(emin, emax, num_points)
    spec_vals = spec(E)
    numerator = xInterpolatedUnivariateSpline(E, spec_vals * pol_func(E))
    denominator = xInterpolatedUnivariateSpline(E, spec_vals)
    return numerator.integral(emin, emax) / denominator.integral(emin, emax)



def broadband_pol_deg(spec, pol_deg, emin=2., emax=8., num_points=500):
    """Calculate the broadband polarization degree for a given spectrum and
    polarization degree vs. energy.
    """
    return _broadband_calc_base(spec, pol_deg, emin, emax, num_points)



def broadband_pol_ang(spec, pol_ang, emin=2., emax=8., num_points=500, degrees=True):
    """Calculate the broadband polarization angle for a given spectrum and
    polarization angle vs. energy.
    """
    pa = _broadband_calc_base(spec, pol_ang, emin, emax, num_points)
    if degrees:
        pa = numpy.degrees(pa)
    return pa



def harmonic_addition(*params):
    """Basic implementation of the harmonic addition theorem, see
    http://mathworld.wolfram.com/HarmonicAdditionTheorem.html

    Args
    ----
    params : 3-element tuple(s)
        The input set of parameters describing the modulation curves that need
        to be added, in the form of an arbitrary number of 3-element tuples
        (F, m, delta) containing the normalization (i.e., the source flux
        associated to the component) the modulation (i.e., the degree of
        polarization of the source) and the phase.

    You should note that, since the modulation curve is defined in terms of
    cos(2 * phi), rather than cos(phi), all the angles in the calculation are
    multiplied by two, and the final arctan is multiplied by 1/2. Keep this in
    mind when you compare the implementation below with the formula at the
    mathworld link.
    """
    F = 0.
    A_square = 0.
    delta_num = 0.
    delta_den = 0.
    for Fi, mi, deltai in params:
        F += Fi
        Ai = Fi * mi
        delta_num += Ai * numpy.sin(2 * deltai)
        delta_den += Ai * numpy.cos(2 * deltai)
        Ai = Fi * mi
        for Fj, mj, deltaj in params:
            Aj = Fj * mj
            A_square += Ai * Aj * numpy.cos(2 * (deltai - deltaj))
    delta = 0.5 * numpy.arctan2(delta_num, delta_den)
    A = numpy.sqrt(A_square)
    m = A / F
    return F, m, delta


def harmonic_component_addition(*components):
    """Harmonic component addition.
    """
    energy = numpy.array([])
    for (spec, pol_deg, pol_ang) in components:
        energy = numpy.union1d(energy, spec.x)
    F = []
    m = []
    delta = []
    for _x in energy:
        params = ((spec(_x), pol_deg(_x), pol_ang(_x)) for \
                  (spec, pol_deg, pol_ang) in components)
        _F, _m, _delta = harmonic_addition(*params)
        F.append(_F)
        m.append(_m)
        delta.append(_delta)
    F = numpy.array(F)
    m = numpy.array(m)
    delta = numpy.array(delta)
    spec = xInterpolatedUnivariateSpline(energy, F)
    pol_deg = xInterpolatedUnivariateSpline(energy, m)
    pol_ang = xInterpolatedUnivariateSpline(energy, delta)
    return spec, pol_deg, pol_ang


def fourier_series_factory(*params):
    """Simple factory class for a generic Fouries series.

    Given an arbitrary number of coefficients of a generic Fuorier expansion,
    this is returning a function that can be evaluated on an arbitrary grid of
    points.

    Note that the Fourier expansion is done in the pulse-phase space, i.e.,
    the output function is expecting an argument in the [0., 1.[ interval
    (peridically repeating itself.)

    Args
    ----
    params : 2-element tuple(s)
        The input set of Fourier coefficients describing the amplitude and
        phase lag of, in the form of 2-element (amplitude, phase) tuples.

    Example
    -------
    >>> import numpy
    >>> from ixpeobssim.srcmodel.polarization import fourier_series_factory
    >>> factory = fourier_series_factory((0.2, 0.5), (0.3, 0.5))
    >>> x = numpy.linspace(0., 1., 100)
    >>> y = factory(x)
    """
    def fourier_series(phase):
        val = 1.
        for i, (a, lag) in enumerate(params):
            val += a * numpy.cos(2 * numpy.pi * (i + 1) * (phase - lag))
        return val
    return fourier_series


def pulse_pol_from_harmonics_spline(mean_flux, mean_pol, *params):
    """Build phase-dependent luminosity and polarization degree, using Fourier
    harmonics composition. This is important for modeling millisecond pulsar
    emission geometry, see:
    Eq.(45) in K. Viironen and J. Poutanen, A&A 426, 985-997 (2004)

    mean_flux must be a scalar in $cm^{-2}s^{-1}keV{-1}$ units;
    mean_pol must be between zero and one.

    For parameters look at Fourier_harmonics function.
    """
    # Shouldn't it be possible to pass the phase from the outside, here?
    phase = numpy.linspace(0., 1., 100)
    vals = fourier_series_factory(*params)(phase)
    pol_deg = mean_pol * vals
    pp_ = mean_flux * vals
    fmt = dict(xlabel='Phase', ylabel='Flux')
    pulse_profile = xInterpolatedUnivariateSpline(phase, pp_, **fmt)
    fmt = dict(xlabel='Phase', ylabel='Polarization Degree')
    polarization_degree = xInterpolatedUnivariateSpline(phase, pol_deg, **fmt)
    return pulse_profile, polarization_degree



class xPolarizationFieldBase:

    """Virtual base class describing a generic, azimuthally simmetric polarization field.

    Arguments
    ---------
    ra0 : float
        The right ascension of the center of the field in decimal degrees.

    dec0 : float
        The declination of the center of the field in decimal degrees.

    radial_profile : float or callable with a signature radial_profile(r, E, t)
        The radial profile for the polarization degree, expressed as a function
        of the distance from the center in decimal degrees, energy and time.
    """

    def __init__(self, ra0, dec0, radial_profile=1.):
        """Constructor.
        """
        self.ra0 = ra0
        self.dec0 = dec0
        if isinstance(radial_profile, numbers.Number):
            self.radial_profile = lambda r, E, t: numpy.full(r.shape, radial_profile)
        else:
            self.radial_profile = radial_profile

    def _delta(self, ra, dec):
        """Return the 2-dimensional offset in ra and dec (or x and y), in decimal degree,
        between a generic point in the sky and the center of the field.

        .. warning:

           Mind this is using the small-angle approximation.
        """
        dx = (ra - self.ra0) * numpy.cos(numpy.radians(self.dec0))
        dy = dec - self.dec0
        return dx, dy

    def _dist_from_center(self, ra, dec):
        """Return the distance from the center of the field in decimal degrees.
        """
        dx, dy = self._delta(ra, dec)
        return numpy.sqrt(dx**2. + dy**2.)

    def polarization_degree(self, E, t, ra, dec):
        """Return the polarization degree as a function of the standard dynamical
        variables.
        """
        r = self._dist_from_center(ra, dec)
        return self.radial_profile(r, E, t)

    def polarization_angle(self, ra, dec):
        """Do nothing method to be reimplemented in derived classes.
        """
        raise NotImplementedError

    def polarization_degree_model(self):
        """Convenience function adapting the signature of the polarization_angle()
        hook to the needs of xpobssim.
        """
        return lambda E, t, ra, dec: self.polarization_degree(E, t, ra, dec)

    def polarization_angle_model(self):
        """Convenience function adapting the signature of the polarization_angle()
        hook to the needs of xpobssim.
        """
        return lambda E, t, ra, dec: self.polarization_angle(ra, dec)



class xRadialPolarizationField(xPolarizationFieldBase):

    def polarization_angle(self, ra, dec):
        """Return the azimuthal angle in the tangential plane, measured in radians
        from the celestial North, at the position (ra, dec) for a purely radial field
        centered at (ra0, dec0).
        """
        dx, dy = self._delta(ra, dec)
        return numpy.arctan2(dx, dy)



class xTangentialPolarizationField(xPolarizationFieldBase):

    def polarization_angle(self, ra, dec):
        """Return the azimuthal angle in the tangential plane, measured in radians
        from the celestial North, at the position (ra, dec) for a purely tangential field
        centered at (ra0, dec0).
        """
        dx, dy = self._delta(ra, dec)
        return numpy.arctan2(dy, -dx)



class xStokesSkyMap:

    """Class representing a Stokes map in sky coordinates.

    This is the basic interface to polarization models for extended sources.
    The basic idea is that we store two bivariate splines for the Q and U
    Stokes parameters, along with the corresponding WCS.

    The class constructor takes raw numpy arrays of Q and U values sampled
    on a regular rectangular grid, but in fact the class is designed in such a
    way that it should rarely be necessary to instantiate a class object from
    the constructor---the model maps will be typically loaded from images stored
    in FITS files, using the convenience class methods. The class supports
    in a transparent way input in three different forms:

    * Q and U;
    * x and y polarization components;
    * polarization degree and angle.

    The class is callable, and returns the Q and U value at given sky
    coordinates (ra and dec).

    Example
    -------
    >>> map_ = xStokesSkyMap.load_from_qu(qpath, upath)
    >>> map_ = xStokesSkyMap.load_from_pda(pdpath, papath)
    """

    def __init__(self, qdata, udata, wcs_, input_file_paths, input_labels):
        """Constructor.
        """
        self.qdata = qdata
        self.udata = udata
        if self.qdata.shape != self.udata.shape:
            logger.info('Shape mismatch: Q is %s, U is %s',
                        self.qdata.shape, self.udata.shape)
            abort('Input arrays to %s must have the same shape' %\
                        self.__class__.__name__)
        self.shape = self.qdata.shape
        self.input_file_paths = input_file_paths
        self.input_labels = input_labels
        self.wcs = wcs_
        # Cache the center and approximate radius of the underlying WCS for
        # future use (e.g., when plotting).
        self.center = self.wcs_center(self.wcs)
        self.radius = self.wcs_radius(self.wcs)
        logger.info('%s', self.wcs)
        # Create the underlying splines.
        x, y = self.map_grid(self.shape)
        self.qspline = xInterpolatedBivariateSplineLinear(x, y, self.qdata)
        self.uspline = xInterpolatedBivariateSplineLinear(x, y, self.udata)
        # Cache the minimum and maximum values of the interpolating grid (in
        # pixel space) as this defines the boundaries of thr grid data for the
        # three-dimensional interpolation.
        self.xmin = x[0]
        self.xmax = x[-1]
        self.ymin = y[0]
        self.ymax = y[-1]

    @staticmethod
    def map_grid(shape):
        """Return grid for map interpolation.

        The spline is performed in logical space, interpolating the values at
        the center of the bins, and the transformation from sky coordinates is
        handled at __call__() time.

        Warning
        -------
        This is tricky, as the pixels in the FITS maps start from 1, while the
        underlying numpy array is zero-indexed. Ideally we want to take the
        interpolation at the center of the bin, and therefore the necessary
        offset is -1 + 0.5 = -0.5. (All of this TBC.)

        See ixpeobssim.srcmodel.img.xFitsImage.__call__() for more information.
        """
        xbins, ybins = shape
        x = numpy.linspace(0, xbins, xbins) - 0.5
        y = numpy.linspace(0, ybins, ybins) - 0.5
        return x, y

    @staticmethod
    def read_map(file_path):
        """Read an input FITS file containg map data.

        This is a generic interface to parse arbitrary bidimensional data
        (depending on the context, these can be in the U, Q or x, y projections
        or polarization degree, angle space).

        The FITS file is supposed to contain an image as the first extension,
        along with the WCS information.
        """
        logger.info('Reading data from %s', file_path)
        wcs_ = wcs.WCS(file_path)
        with fits.open(file_path) as hdu_list:
            # Mind we transpose the image to fit into a numpy array with the
            # proper orientation.
            data = hdu_list[0].data.transpose()
        return data, wcs_

    @staticmethod
    def pixel_shape(wcs_):
        """Small convenience function to retrieve the pixel shape of a WCS object.

        We need this because the private attributes "_naxis1" and "_naxis2" have
        been deprecated since astropy version 3.1 in favor of the  "pixel_shape"
        property.
        """
        if ASTROPY_VERSION <= '3.1':
            return wcs_._naxis1, wcs_._naxis2
        return wcs_.pixel_shape

    @staticmethod
    def compare_wcs(wcs1, wcs2):
        """Compare two WCS objects to make sure they are identical.

        I am sure we can do this in a better way, but the way the information
        is stored in the WCS seems really convoluted. At the zero order we are
        trying to make sure that two WCS objects have the same center, axes and
        extension.
        """
        # Test some "simple" properties...
        for prop in ['naxis', 'crpix', 'crval']:
            val1 = getattr(wcs1.wcs, prop)
            val2 = getattr(wcs2.wcs, prop)
            diff = val1 - val2
            if isinstance(diff, numpy.ndarray):
                diff = diff.all()
            if diff:
                abort('WCS mismatch for %s (%s vs %s)' % (prop, val1, val2))
        # ...then the shape...
        shape1 = xStokesSkyMap.pixel_shape(wcs1)
        shape2 = xStokesSkyMap.pixel_shape(wcs2)
        if shape1[0] != shape2[0] or shape1[1] != shape2[1]:
            abort('WCS mismatch for pixel_shape (%s vs %s)' % (shape1, shape2))
        # and finally the delta.
        if wcs1.wcs.has_cd():
            delta1 = abs(wcs1.wcs.cd.diagonal())
            delta2 = abs(wcs2.wcs.cd.diagonal())
        else:
            delta1 = abs(wcs1.wcs.cdelt)
            delta2 = abs(wcs2.wcs.cdelt)
        diff = (delta1 - delta2).all()
        if diff:
            abort('WCS mismatch for delta (%s vs %s)' % (delta1, delta2))

    @staticmethod
    def wcs_center(wcs_):
        """Return the coordinates of the center of a given WSC object.
        """
        ra, dec = [float(val) for val in wcs_.wcs.crval]
        return ra, dec

    @staticmethod
    def wcs_radius(wcs_, default=arcmin_to_degrees(6.)):
        """Return the approximate radius of a given WCS object.

        The FITS specification dictates that when CD is present, CDELT should
        be set to 1, and is effectively ignored, see
        https://docs.astropy.org/en/stable/api/astropy.wcs.Wcsprm.html
        """
        pixel_shape = xStokesSkyMap.pixel_shape(wcs_)
        # If the pixel shape is None there is nothing to do other than returning
        # a default value comparable with the field of view of the instrument.
        if pixel_shape is None:
            return default
        # Retrieve the pixel size based on the appropriate method for the WCS
        # at hand.
        if wcs_.wcs.has_cd():
            delta = abs(wcs_.wcs.cd.diagonal())
        else:
            delta = abs(wcs_.wcs.cdelt)
        # Calculate the radius. It is not entirely clear to me why in this
        # context we need to divide for the cosine of the declinations in
        # both Ra and Dec, but this seem to yield sensible results.
        ra, dec = xStokesSkyMap.wcs_center(wcs_)
        radius = 0.5 * delta * pixel_shape / numpy.cos(numpy.radians(dec))
        return radius.max()

    @classmethod
    def load_from_qu(cls, q_file_path, u_file_path):
        """Load a sky map of Stokes parameters from Q and U data.
        """
        qdata, qwcs = cls.read_map(q_file_path)
        udata, uwcs = cls.read_map(u_file_path)
        cls.compare_wcs(qwcs, uwcs)
        labels = ('Q', 'U')
        return cls(qdata, udata, qwcs, (q_file_path, u_file_path), labels)

    @classmethod
    def load_from_pda(cls, pd_file_path, pa_file_path):
        """Load a sky map of Stokes parameters from polarization degree and
        angle.
        """
        pddata, pdwcs = cls.read_map(pd_file_path)
        padata, pawcs = cls.read_map(pa_file_path)
        cls.compare_wcs(pdwcs, pawcs)
        qdata = xModelStokesParameters.q(pddata, padata)
        udata = xModelStokesParameters.u(pddata, padata)
        labels = ('Polarization degree', 'Polarization angle')
        return cls(qdata, udata, pdwcs, (pd_file_path, pa_file_path), labels)

    def __call__(self, ra, dec):
        """Overloaded method.

        This interpolates the Stokes parameters in sky coordinates.

        Warning
        -------
        At this point we don't have any sensible mechanism to extrapolate when
        the input coordinates happen to be outside the underlying interpolating
        grid, and we resort set the Stokes parameters to zero in that case,
        issuing a warning message.
        """
        # Convert to pixel coordinates for interpolation.
        x, y = self.wcs.wcs_world2pix(ra, dec, 0)
        # Flag the points outside the interpolation grid.
        _xmask = numpy.logical_or(x < self.xmin, x > self.xmax)
        _ymask = numpy.logical_or(y < self.ymin, y > self.ymax)
        pos_mask = numpy.logical_or(_xmask, _ymask)
        q = self.qspline(x, y)
        u = self.uspline(x, y)
        # Explicitely set to zero the Stokes parameters for the points outside
        # the interpolator domain.
        q[pos_mask] = 0.
        u[pos_mask] = 0.
        return q, u

    def polarization_vector(self, ra, dec):
        """Return the x and y components of the polarization vector at given
        sky coordinates. .
        """
        return xModelStokesParameters.qu_to_xy(*self(ra, dec))

    def polarization_degree(self, ra, dec):
        """Return the polarization degree value given a ra and dec.
        """
        q, u = self(ra, dec)
        return xModelStokesParameters.polarization_degree(q, u)

    def polarization_angle(self, ra, dec):
        """Return the polarization angle value given a ra and dec.
        """
        q, u = self(ra, dec)
        return xModelStokesParameters.polarization_angle(q, u)

    def polarization_degree_model(self):
        """Convenience method to adapt the signature of the __call__ class
        method to the one needed at simulation time.
        """
        return lambda E, t, ra, dec: self.polarization_degree(ra, dec)

    def polarization_angle_model(self):
        """Convenience method to adapt the signature of the __call__ class
        method to the one needed at simulation time.
        """
        return lambda E, t, ra, dec: self.polarization_angle(ra, dec)

    def plot_input_data(self, suffix=None):
        """Plot the underlying input arrays used to build the sky map.
        """
        figs = []
        for file_path, label in zip(self.input_file_paths, self.input_labels):
            figure_name = label
            if suffix is not None:
                figure_name = '%s - %s' % (label, suffix)
            plt.figure(figure_name)
            img = xFITSImage(file_path)
            figs.append(img.plot(zlabel=label))
        return figs

    def _make_plot(self, data, label, vmin, vmax, cmap):
        """Convenience function to factor out the code in common to the
        different plotting routines.

        Mind we have to transpose the underlying array for the plot to work.
        """
        data = data.transpose()
        plt.subplot(projection=self.wcs)
        plt.imshow(data, origin='lower', vmin=vmin, vmax=vmax, cmap=cmap)
        plt.grid(color='white', ls='dashed')
        plt.xlabel('RA (J2000)')
        plt.ylabel('Dec (J200)')
        plt.colorbar(label=label)

    def plot_polarization_degree(self, vmin=None, vmax=None, cmap='gist_heat'):
        """Plot the polarization degree as a function of the sky direction.
        """
        data = xModelStokesParameters.polarization_degree(self.qdata, self.udata)
        label = 'Polarization degree'
        self._make_plot(data, label, vmin, vmax, cmap)

    def plot_polarization_angle(self, vmin=None, vmax=None, cmap='gist_heat', degrees=True):
        """Plot the polarization angle as a function of the sky direction.
        """
        data = xModelStokesParameters.polarization_angle(self.qdata, self.udata)
        if degrees is True:
            data = numpy.degrees(data)
            units = 'degrees'
        else:
            units = 'rad'
        label = 'Polarization angle [%s]' % units
        self._make_plot(data, label, vmin, vmax, cmap)

    def plot_q(self, vmin=None, vmax=None, cmap='gist_heat'):
        """Plot Q as a function of the sky direction.
        """
        self._make_plot(self.qdata, 'Q', vmin, vmax, cmap)

    def plot_u(self, vmin=None, vmax=None, cmap='gist_heat'):
        """Plot U as a function of the sky direction.
        """
        self._make_plot(self.udata, 'Q', vmin, vmax, cmap)

    def plot_arrows(self, nside=25, ra0=None, dec0=None, radius=None,
                    threshold=0.0, **kwargs):
        """Overlay the polarization arrows on top of a map.
        """
        if ra0 is None:
            ra0 = self.center[0]
        if dec0 is None:
            dec0 = self.center[1]
        if radius is None:
            radius = self.radius
        else:
            radius = arcmin_to_degrees(radius)
        grid = square_sky_grid(nside, (ra0, dec0), radius)
        plot_arrows(grid, self.polarization_vector, threshold, **kwargs)



class xStokesSkyCubeLayer(xStokesSkyMap):

    """Small convenience class representing a layer of an xStokesSkyCube object.

    This is essentially a lightweight wrapper of the xStokesSkyMap class, on top
    of which we add a few extra members (e.g., the energy bounds) that are
    necessary when a map is used in the context of a cube.
    """

    # pylint: disable=attribute-defined-outside-init

    def set_energy_range(self, emin, emax=None):
        """Set the minimum and maximum energy for the layer.

        If emax is None it is assumed that the layer has been calculated at
        a specified energy (as opposed to an energy range), and is treated as
        such---the minimum and maximum energies are the same.
        """
        if emax is None:
            emax = emin
        self.emin = emin
        self.emax = emax
        self.emean = 0.5 * (self.emin + self.emax)

    def energy_range_label(self):
        """Return a text label encoding the proper energy range for the layer.
        """
        if self.emin == self.emax:
            return '%.2f keV' % self.emin
        return '%.2f--%.2f keV' % (self.emin, self.emax)

    def label(self):
        """Return a generic text label for the layer.
        """
        return 'Stokes sky-cube layer (%s)' % self.energy_range_label()



class xStokesSkyCube:

    """Class describing a Stoke sky-cube, i.e., a coherent collection of Stokes
    sky maps in different energy layers.
    """

    def __init__(self):
        """Constructor.
        """
        self.__layers = []
        self.wcs = None
        self.shape = None
        self.center = None
        self.__qinterpolator = None
        self.__uinterpolator = None

    def layer(self, i):
        """Return a generic cube layer by index.
        """
        return self.__layers[i]

    def xy_grid_bounds(self):
        """Return the spatial bounds (in pixel space) for the underlying
        interpolating grid.

        Note that, since we make sure at fill time that all the layers have the
        same WCS, here we take the bounds from the first layer.
        """
        layer = self.__layers[0]
        return layer.xmin, layer.xmax, layer.ymin, layer.ymax

    def _add_layer(self, layer, emin, emax=None):
        """Add a layer to the cube.

        This is a generic function that is factoring all the operations that
        are performed *after* a layer is created (e.g., updating the WCS
        information and making sure that the WCS and underlying data are
        consistent across layers). Users should typically *not* call this
        method (which is why we make it provate), but use one of the methods
        below, that create and add a new layer from suitable FITS files.
        """
        # If this is the first layer we are adding, take note of the WCS and
        # the shape of the underlying maps...
        if self.wcs is None:
            self.wcs = layer.wcs
            self.center = xStokesSkyMap.wcs_center(self.wcs)
            self.shape = layer.shape
        # ...otherwise make sure that all follwing layers have a WCS and a shape
        # which are consistent with the first one.
        else:
            xStokesSkyMap.compare_wcs(self.wcs, layer.wcs)
            assert layer.shape == self.shape
        layer.set_energy_range(emin, emax)
        self.__layers.append(layer)

    def add_layer_qu(self, q_file_path, u_file_path, emin, emax=None):
        """Add a layer for Q and U maps.
        """
        layer = xStokesSkyCubeLayer.load_from_qu(q_file_path, u_file_path)
        self._add_layer(layer, emin, emax)

    def add_layer_pda(self, pd_file_path, pa_file_path, emin, emax=None):
        """Add a layer for polarization degree and angle.
        """
        layer = xStokesSkyCubeLayer.load_from_pda(pd_file_path, pa_file_path)
        self._add_layer(layer, emin, emax)

    def _calculate_grid_data(self):
        """Calculate all the data structures needed to interpolate the
        polarization metrics.
        """
        # Prepare the binning for the two interpolator objects.
        x, y = xStokesSkyMap.map_grid(self.shape)
        z = numpy.array([layer.emean for layer in self.__layers])
        # Fill the data from the underlying layers.
        qdata = numpy.zeros((*self.shape, *z.shape))
        udata = numpy.zeros((*self.shape, *z.shape))
        for k, layer in enumerate(self.__layers):
            qdata[:, :, k] = layer.qdata
            udata[:, :, k] = layer.udata
        # Create the actual interpolators.
        self.__qinterpolator = RegularGridInterpolator((x, y, z), qdata)
        self.__uinterpolator = RegularGridInterpolator((x, y, z), udata)

    def __call__(self, ra, dec, energy):
        """Return the interpolated Stokes parameters at given sky coordinates
        and energies.

        Warning
        -------
        At this point we don't have any sensible mechanism to extrapolate when
        the input coordinates happen to be outside the underlying interpolating
        grid, and we resort set the Stokes parameters to zero in that case,
        issuing a warning message.
        """
        if self.__qinterpolator is None or self.__uinterpolator is None:
            self._calculate_grid_data()
        x, y = self.wcs.wcs_world2pix(ra, dec, 0)
        # Horrible hack to avoid an interpolation error.
        # Deal with sky coordinates, first...
        xmin, xmax, ymin, ymax = self.xy_grid_bounds()
        _xmask = numpy.logical_or(x < xmin, x > xmax)
        _ymask = numpy.logical_or(y < ymin, y > ymax)
        pos_mask = numpy.logical_or(_xmask, _ymask)
        n = pos_mask.sum()
        if n > 0:
            frac = 100. * n / float(len(x))
            logger.error(INTERPOLATION_POS_BOUND_ERROR_MSG, n, frac)
            logger.info('Full input Ra range: %.5f--%.5f deg', ra.min(), ra.max())
            logger.info('Full input Dec range: %.5f--%.5f deg', dec.min(), dec.max())
            x = numpy.clip(x, xmin, xmax)
            y = numpy.clip(y, ymin, ymax)
        # ... and then with energies.
        emin = self.__layers[0].emin
        emax = self.__layers[-1].emax
        ene_mask = numpy.logical_or(energy < emin, energy > emax)
        n = ene_mask.sum()
        if n > 0:
            frac = 100. * n / float(len(energy))
            logger.error(INTERPOLATION_ENE_BOUND_ERROR_MSG, n, frac)
            logger.info('Full input energy range: %.3f--%.3f keV',
                        energy.min(), energy.max())
            energy = numpy.clip(energy, emin, emax)
        # End of the hack.
        q = self.__qinterpolator((x, y, energy))
        u = self.__uinterpolator((x, y, energy))
        # Set to zero both Stokes parameters for all the positions ending
        # up outside the domain of the underlying interpolator.
        q[pos_mask] = 0.
        u[pos_mask] = 0.
        return q, u

    def polarization_degree(self, ra, dec, energy):
        """Return the polarization degree value given a ra, dec and energy.
        """
        q, u = self(ra, dec, energy)
        return xModelStokesParameters.polarization_degree(q, u)

    def polarization_angle(self, ra, dec, energy):
        """Return the polarization angle value given a ra, dec and energy.
        """
        q, u = self(ra, dec, energy)
        return xModelStokesParameters.polarization_angle(q, u)

    def polarization_degree_model(self):
        """Convenience method to adapt the signature of the __call__ class
        method to the one needed at simulation time.
        """
        return lambda E, t, ra, dec: self.polarization_degree(ra, dec, E)

    def polarization_angle_model(self):
        """Convenience method to adapt the signature of the __call__ class
        method to the one needed at simulation time.
        """
        return lambda E, t, ra, dec: self.polarization_angle(ra, dec, E)

    def plot(self, **kwargs):
        """Plot all the layers.
        """
        for layer in self.__layers:
            plt.figure('Polarization degree - %s' % layer.label())
            layer.plot_polarization_degree(**kwargs)
            plt.figure('Polarization angle - %s' % layer.label())
            layer.plot_polarization_angle(**kwargs)

    def plot_input_data(self):
        """Plot the input data for all the layers.
        """
        for i, layer in enumerate(self.__layers):
            print('layer %d (%s)' % (i, layer.energy_range_label()))
            layer.plot_input_data(suffix=layer.label())

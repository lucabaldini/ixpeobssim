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

"""Model component and region of interest definition.
"""

from __future__ import print_function, division

from collections import OrderedDict
import numbers

import numpy
import scipy.stats
from astropy.io import fits
from astropy import wcs

from ixpeobssim.srcmodel.img import xFITSImage
from ixpeobssim.srcmodel.spectrum import xCountSpectrum
from ixpeobssim.srcmodel.ephemeris import phase_function, time_list
from ixpeobssim.evt.event import xEventList
from ixpeobssim.evt.fmt import standard_radec_to_xy
from ixpeobssim.evt.ixpesim import xPhotonList, build_tow_response
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.core.stokes import xModelStokesParameters
from ixpeobssim.utils.units_ import keV_to_erg, ergcms_to_mcrab, degrees_to_arcmin, erg_to_keV
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.astro import ds9_region_filter_sky, angular_separation
from ixpeobssim.utils.astro import wcs_digitize, wcs_to_sky_meshgrid
from ixpeobssim.utils.time_ import met_to_mjd
from ixpeobssim.instrument.gpd import phi_to_detphi
from ixpeobssim.instrument.mma import sky_to_gpd, apply_dithering, parse_dithering_kwargs
from ixpeobssim.instrument import gpd
# pylint: disable=consider-using-from-import
import ixpeobssim.utils.chandra as chandra

# pylint: disable=invalid-name, too-many-arguments, too-many-locals, too-many-lines
# pylint: disable=too-many-instance-attributes, no-member


class xModelComponentBase:

    """Base class for all the model components (i.e., the sources) populating
    a generic region of interest for a simulation.

    Historically this base class was tailored to celestial sources, but since
    different use cases emerged along the way that do not really fit in the
    original picture (e.g., instrumental background and calibration flat fields,
    where everything happens within the detector, and there is no concept of
    mirror effective area nor vignetting, just to name a few) the decision was
    taken to split the very fundamental features of the object into this minimal
    base class and have a separate base class for the celestial objects.

    In addition to some fundamental members and method, the class includes a
    handful of static methods that encapsulate useful operations of general
    interest, such as generating event times and energy according to simple
    analytical distributions.

    Arguments
    ---------
    name : string
        The name of the source.

    identifier : int, optional
        A unique identifier of the source within a ROI model. Defaults to None
        and is, generally, automatically assigned when building the ROI (i.e.,
        you don't have to worry about it, but it's handy to have in the
        output event list).
    """

    def __init__(self, name, identifier=None):
        """Constructor.
        """
        self.name = name
        self.identifier = identifier

    def rvs_event_list(self, parent_roi, irf_set, **kwargs):
        """Return an event list for the model component.

        Arguments
        ---------
        parent_roi : xROIModel instance
            The parent region of interest (i.e., the xROIModel instance the
            model component belongs to).

        irf_set : xIRFSet instance
            The set of detector response functions to be used for the
            convolution of the physical quantities.

        **kargs : keyword arguments
            At runtime this is precisely the complete set of command-line
            options that xpobssim is run with (i.e., you do get all the goodies
            describing the configuration of the simulation, such as the start
            mission elapsed time and the simulation duration).

        Note
        ----
        This is a do-nothing function that should be re-implemented in every
        derived class.

        This is also the main xModelComponentBase method, since all the Physics
        of the model itself is encapsulated here. In addition, the method is
        aware of the parent region of interest object and the instrument
        response functions, which are essential ingredients to create the output
        event list.

        As a rule of thumb:

        * all the operations that makes sense to perform on a component by
          component basis (e.g., the convolution with the instrument response
          functionsm, including the vignetting) happen in the body of this
          function;
        * all the operations (such as the application of the deadtime) that need
          to be performed on the final event list after all the model model
          components have been merged togethere are deferred to the downstream
          code.
        """
        raise NotImplementedError

    def rvs_photon_list(self, parent_roi, irf_set, **kwargs):
        """Return a photon list for the model component.
        """
        raise NotImplementedError

    @staticmethod
    def poisson(average):
        """Return a number extracted from a Poisson distribution with a given
        average.
        """
        return numpy.random.poisson(average)

    @staticmethod
    def uniform_time(size, start_met, duration):
        """Return a sorted array of a given size of event times uniformly
        distributed between a minimum and a maximum time.
        """
        time_ = numpy.random.uniform(start_met, start_met + duration, size)
        time_.sort()
        return time_

    @staticmethod
    def uniform_phi(size):
        """Return an array of a given size of azimuthal angles distributed
        between -pi and pi.

        This can be used to replace the fully-fledged xAzimuthalResponse
        generator in cases where one knows since the beginning that a flat
        distribution is needed (it goes without saying that this is much less
        demanding from a computational standpoint).
        """
        return numpy.random.uniform(-numpy.pi, numpy.pi, size)

    @staticmethod
    def uniform_rectangle(size, half_side_x, half_side_y):
        """Return two arrays of a given size with random x and y coordinates
        uniformly distributed within a rectangle of a given half-sides centered
        in 0, 0.
        """
        detx = numpy.random.uniform(-half_side_x, half_side_x, size)
        dety = numpy.random.uniform(-half_side_y, half_side_y, size)
        return detx, dety

    @staticmethod
    def uniform_square(size, half_side):
        """Specialized function generating random coordinates uniformly
        distributed within a square.
        """
        return self.uniform_square(size, half_side, half_side)

    @classmethod
    def convolve_sky_direction(cls, mc_ra, mc_dec, parent_roi, psf):
        """Smear the sky position with the PSF and calculate the corresponding
        digitized x and y values.
        """
        ra, dec = psf.smear(mc_ra, mc_dec)
        x, y = standard_radec_to_xy(ra, dec, parent_roi.ra, parent_roi.dec)
        return ra, dec, x, y

    @staticmethod
    def parse_time_kwargs(**kwargs):
        """Parse the keyword arguments related to the duration of the observation.
        """
        return kwargs.get('start_met'), kwargs.get('duration')

    def __str__(self):
        """String formatting.
        """
        return '{} "{}" (id = {})'.format(self.__class__.__name__, self.name,
                                          self.identifier)



class xCelestialModelComponentBase(xModelComponentBase):

    """Base class for the source object.

    Note that the source identifier defaults to none and is typically assign
    after the fact when the source itself is added to a source model.

    Arguments
    ---------
    name : string
        The name of the source.

    photon_spectrum : function
        The function object representing the photon spectrum.

    polarization_degree : function
        The function object representing the polarization degree.

    polarization_angle : function
        The function object representing the polarization angle.

    column_density : float
        The value of the column density (in cm^-2) used to calculate the
        Galactic absorption. Defaults to 0.

    redshift : float
        The source redshift. Defaults to 0.

    identifier : int
        A unique identifier of the source within a ROI model. Defaults to None
        and is, generally, automatically assigned when building the ROI (i.e.,
        you don't have to worry about it, but it's handy to have in the
        output event list).
    """

    def __init__(self, name, photon_spectrum, polarization_degree,
                 polarization_angle, column_density=0., redshift=0.,
                 identifier=None):
        """Constructor.
        """
        xModelComponentBase.__init__(self, name, identifier)
        self.setup(photon_spectrum, polarization_degree, polarization_angle)
        self.column_density = column_density
        self.redshift = redshift
        # Now the advanced settings for the sampling time grid and the order of
        # the bivariate spline for the count spectrum. See
        # https://bitbucket.org/ixpesw/ixpeobssim/issues/55
        # for more discussion about this.
        self._count_spectrum_ny = 200
        self._count_spectrum_kx = 3
        self._count_spectrum_ky = 3
        # We cache here the integral flux, since we don't want to (potentially)
        # recalculate the integral each time the __str__() method is called.
        emin = 2.
        emax = 8.
        t = 0.
        self.flux_ergcms = self.calculate_integral_flux(emin, emax, t)
        self.flux_mcrab = ergcms_to_mcrab(self.flux_ergcms)
        self.flux_label = 'Unabsorbed flux @ t = %d: %.3e erg/cm2/s' % (t, self.flux_ergcms)
        self.flux_label += ' (%.2f mcrab)' % self.flux_mcrab

    def setup(self, photon_spectrum, polarization_degree, polarization_angle):
        """Setup the model component in terms of photon spectrum and
        polarization degree and angle.
        """
        self.photon_spectrum = photon_spectrum
        self.polarization_degree = polarization_degree
        self.polarization_angle = polarization_angle

    def set_count_spectrum_params(self, ny=200, kx=3, ky=3):
        """Adjust the parameters used at runtime for the creation of the
        count spectrum.

        This is the main hook to control how the count spectrum for a given
        source component is created, see
        https://bitbucket.org/ixpesw/ixpeobssim/issues/55
        """
        self._count_spectrum_ny = ny
        self._count_spectrum_kx = kx
        self._count_spectrum_ky = ky

    def calculate_integral_flux(self, emin=2., emax=8., t=0., num_points=250, erg=True):
        """Return the integral source flux at a generic time.

        This is achieved by taking a "slice" of the source spectrum
        at that time and integrating between a minimum and maximum energy.

        Arguments
        ---------
        emin : float
            The minimum integration energy (default 2 keV).

        emax : float
            The maximum integration energy (default 8 keV).

        t : float
            The time (default is 0).
        """
        x = numpy.linspace(emin, emax, num_points)
        y = x * self.photon_spectrum(x, t)
        flux = xInterpolatedUnivariateSplineLinear(x, y).integral(emin, emax)
        if erg:
            flux = keV_to_erg(flux)
        return flux

    def calculate_average_polarization(self, egrid, tgrid, degrees=False):
        """Calculate the average polarization degree and angle in a given
        energy and time (or phase) interval.

        Warning
        -------
        This is deprecated in favor of the model_pol_average facilities.
        """
        # If either of the inputs is a number, we need to turn it into the
        #
        if isinstance(egrid, numbers.Number):
            egrid = numpy.array([egrid], float)
        if isinstance(tgrid, numbers.Number):
            tgrid = numpy.array([tgrid], float)
        energy, time_ = numpy.meshgrid(egrid, tgrid)
        spec = self.photon_spectrum(energy, time_)
        pol_deg = self.polarization_degree(energy, time_, None, None)
        pol_ang = self.polarization_angle(energy, time_, None, None)
        q = xModelStokesParameters.q(pol_deg, pol_ang)
        u = xModelStokesParameters.u(pol_deg, pol_ang)
        norm = spec.sum()
        Q = (q * spec).sum() / norm
        U = (u * spec).sum() / norm
        pol_deg = xModelStokesParameters.polarization_degree(Q, U)
        pol_ang = xModelStokesParameters.polarization_angle(Q, U)
        if degrees:
            pol_ang = numpy.degrees(pol_ang)
        return pol_deg, pol_ang

    def build_intensity_map(self, wcs_, num_points=1000000):
        """Base function to calculate the underlying intensity map over a regular
        sky grid.

        In the base class this is implemented via brute force, i.e., we throw a
        bunch or random sky directions via a rvs_sky_coordinates() call and then
        we bin them on the approriate wcs object. In reality, for most of the
        source classes, this can be done much more effectively by other means,
        with two major advantages:

        * the calculation can be made faster;
        * the numerical noise due to the Monte Carlo integration can be entirely
          avoided.

        Note that, whenever possible, this function is supposed to return the
        actual model values on a grid of points, rather than a brute-force
        Monte Carlo approximation. In practice, we should strive to reimplement
        this method in all sub-classes, with the only exception of extended
        sources based on FITS images.
        """
        ra, dec = self.rvs_sky_coordinates(num_points)
        data = wcs_digitize(wcs_, ra, dec)
        data /= data.sum()
        return data

    def __str__(self):
        """String formatting.

        Note that this base class has no ra and dec attributes, since they're
        not really relevant for extended sources initiazed from FITS images,
        but since all the other derived classed do have such attributes,
        check whether that's the case and, if necessary, print out the
        source position. This avoids code duplication in the derived classes.
        """
        text = xModelComponentBase.__str__(self)
        text += '\n    Galactic column density: %.3e cm^{-2}' %\
                self.column_density
        text += '\n    Redshift: %.3f' % self.redshift
        text += '\n    %s' % self.flux_label
        if hasattr(self, 'ra') and hasattr(self, 'dec'):
            text += '\n    Position: RA = %s deg, Dec = %s deg' %\
                    (self.ra, self.dec)
        return text

    def rvs_sky_coordinates(self, size=1):
        """Generate random coordinates for the model component.

        This is a do-nothing function and should be re-implemented by
        each derived class.

        Arguments
        ---------
        size : float
            The number of sky coordinate pairs to be generated.
        """
        raise NotImplementedError

    def sampling_time_grid(self, start_met, duration):
        """Return the time grid used to sample the lightcurve when generating
        events.
        """
        return numpy.linspace(start_met, start_met + duration, self._count_spectrum_ny)

    def create_count_spectrum(self, aeff, time_grid, **kwargs):
        """Create the count spectrum object at the base of the source
        simulation.

        This is the bit where we convolve the source spectrum with the effective
        area.

        Args
        ----
        aeff : xEffectiveArea object
            The effective area as a function of the energy---this can be a generic
            callable accepting the energy as the first (and only) argument.

        time_grid : array_like
            The time (or phase) grid over which the count spectrum is created.
        """
        _kwargs = dict(emin=kwargs.get('emin'), emax=kwargs.get('emax'),
            kx=self._count_spectrum_kx, ky=self._count_spectrum_ky)
        args = self.photon_spectrum, aeff, time_grid, self.column_density, self.redshift
        return xCountSpectrum(*args, **_kwargs)

    def _rvs_phi(self, modf, energy, time_, ra, dec):
        """Generate random photoelectron directions with the proper
        distribution.

        This is evaluating the polarization degree and angle at the
        appropriate energy, time and position in the sky and returning a
        corresponding array of photoelectron directions.
        """
        pd = self.polarization_degree(energy, time_, ra, dec)
        pa = self.polarization_angle(energy, time_, ra, dec)
        return modf.rvs_phi(energy, pd, pa)

    def _rvs_seed_event_list(self, parent_roi, irf_set, **kwargs):
        """Generate an event list with random time and energy values for the
        model component and for a given effective area and simulation setup.
        """
        # Create the count spectrum.
        time_grid = self.sampling_time_grid(*self.parse_time_kwargs(**kwargs))
        count_spectrum = self.create_count_spectrum(irf_set.aeff, time_grid, **kwargs)
        # Extract the number of events to be generated.
        num_events = self.poisson(count_spectrum.light_curve.norm())
        logger.info('About to generate %d events...', num_events)
        # Extract the event times.
        time_ = count_spectrum.rvs_event_times(num_events)
        # Apply the GTIs.
        time_, _ = kwargs.get('gti_list').filter_event_times(time_)
        num_events = len(time_)
        if num_events == 0:
            return xEventList()
        event_list = xEventList(time_, self.identifier)
        # Extract the event energies and perform the pha analysis.
        mc_energy = count_spectrum.rvs(time_)
        mc_pha, mc_pi = irf_set.edisp.pha_analysis(mc_energy)
        # Extract the sky coordinates.
        mc_ra, mc_dec = self.rvs_sky_coordinates(num_events)
        mc_x, mc_y = standard_radec_to_xy(mc_ra, mc_dec, parent_roi.ra, parent_roi.dec)
        # Extract the photoelectron emission directions.
        phi = self._rvs_phi(irf_set.modf, mc_energy, time_, mc_ra, mc_dec)
        # Rotate the phi angles in the gpd reference frame
        detphi = phi_to_detphi(phi, irf_set.du_id, kwargs.get('roll'))
        cols = [mc_energy, mc_pha, mc_pi, mc_ra, mc_dec, mc_x, mc_y, phi, detphi]
        event_list.set_seed_columns(*cols)
        return event_list

    @classmethod
    def convolve_event_list(cls, event_list, parent_roi, irf_set, **kwargs):
        """Convolve a Monte Carlo event list with a set of instrument
        response functions.

        Note that we factor this functionality as a separate class method
        so that it can be used in different context, e.g., in the Chandra to
        IXPE converter.
        """
        # If the input event list is empty, just forward it out.
        if event_list.empty():
            return event_list
        # Retrieve the Monte Carlo truth from the original event list.
        mc_energy = event_list.mc_energy()
        mc_ra, mc_dec = event_list.mc_sky_coordinates()
        # Convolve the energy with the energy dispersion.
        energy, pha, pi = irf_set.edisp.convolve_energy(mc_energy)
        # Smear the sky position with the PSF.
        ra, dec, x, y = cls.convolve_sky_direction(mc_ra, mc_dec, parent_roi, irf_set.psf)
        time_ = event_list.time()
        # Apply the dithering effect to the pointing direction (if needed).
        dither_params = parse_dithering_kwargs(**kwargs)
        ra_pnt, dec_pnt = apply_dithering(time_, parent_roi.ra, parent_roi.dec, dither_params)
        # Project the sky positions onto the gpd reference frame.
        detx, dety = sky_to_gpd(ra, dec, time_, ra_pnt, dec_pnt, irf_set.du_id, kwargs.get('roll'))
        # Fill the actual event list.
        cols = [pha, pi, energy, ra, dec, x, y, detx, dety]
        event_list.set_rec_columns(*cols)
        # If the vignetting is enabled, apply it.
        if kwargs.get('vignetting'):
            event_list.apply_vignetting(irf_set.vign, ra_pnt, dec_pnt)
        return event_list

    def rvs_event_list(self, parent_roi, irf_set, **kwargs):
        """Extract a random event list for the model component.
        """
        event_list = self._rvs_seed_event_list(parent_roi, irf_set, **kwargs)
        self.convolve_event_list(event_list, parent_roi, irf_set, **kwargs)
        return event_list

    def rvs_photon_list(self, parent_roi, irf_set, **kwargs):
        """Extract a random photon list.
        """
        # Build the custom response at the top of the window.
        aeff_spline, _ = build_tow_response(irf_set)
        # Create the count spectrum.
        time_grid = self.sampling_time_grid(*self.parse_time_kwargs(**kwargs))
        count_spectrum = self.create_count_spectrum(aeff_spline, time_grid, **kwargs)
        # Extract the number of events to be generated.
        num_events = self.poisson(count_spectrum.light_curve.norm())
        logger.info('About to generate %d photons...', num_events)
        # Extract the event times.
        time_ = count_spectrum.rvs_event_times(num_events)
        # Apply the GTIs.
        time_, _ = kwargs.get('gti_list').filter_event_times(time_)
        num_events = len(time_)
        if num_events == 0:
            return xPhotonList()
        photon_list = xPhotonList(time_, self.identifier)
        energy = count_spectrum.rvs(time_)
        mc_ra, mc_dec = self.rvs_sky_coordinates(num_events)
        # Smear the sky position with the PSF.
        ra, dec, _, _ = self.convolve_sky_direction(mc_ra, mc_dec, parent_roi, irf_set.psf)
        # Apply the dithering effect to the pointing direction (if needed).
        dither_params = parse_dithering_kwargs(**kwargs)
        roll_angle = kwargs.get('roll')
        ra_pnt, dec_pnt = apply_dithering(time_, parent_roi.ra, parent_roi.dec, dither_params)
        # Project the sky positions onto the gpd reference frame.
        detx, dety = sky_to_gpd(ra, dec, time_, ra_pnt, dec_pnt, irf_set.du_id, roll_angle)
        pol_deg = self.polarization_degree(energy, time_, mc_ra, mc_dec)
        pol_ang = self.polarization_angle(energy, time_, mc_ra, mc_dec)
        # Rotate the polarization angle from the sky reference frame to the
        # GPD reference frame.
        pol_ang = phi_to_detphi(pol_ang, irf_set.du_id, roll_angle)
        photon_list.fill(energy, ra, dec, detx, dety, pol_deg, pol_ang)
        # If the vignetting is enabled, apply it.
        if kwargs.get('vignetting'):
            photon_list.apply_vignetting(irf_set.vign, ra_pnt, dec_pnt)
        return photon_list



class xPointSource(xCelestialModelComponentBase):

    """Class representing a steady point source.

    See :py:class:`ixpeobssim.srcmodel.roi.xCelestialModelComponentBase` for the
    signature of the base class.

    Arguments
    ---------
    ra : float
        The right ascension of the source (in decimal degrees).

    dec : float
        The declination of the source (in decimal degrees).
    """

    def __init__(self, name, ra, dec, photon_spectrum, polarization_degree,
                 polarization_angle, column_density=0., redshift=0.,
                 identifier=None):
        """Constructor.
        """
        args = name, photon_spectrum, polarization_degree, polarization_angle,\
            column_density, redshift, identifier
        xCelestialModelComponentBase.__init__(self, *args)
        self.ra = ra
        self.dec = dec

    # pylint: disable=arguments-differ
    def build_intensity_map(self, wcs_):
        """Overloaded method.
        """
        return wcs_digitize(wcs_, self.ra, self.dec)

    def rvs_sky_coordinates(self, size=1):
        """Generate random coordinates for the model component.

        This is returning an array of the proper length with identical values.

        Arguments
        ---------
        size : float
            The number of sky coordinate pairs to be generated.
        """
        ra = numpy.full(size, self.ra)
        dec = numpy.full(size, self.dec)
        return (ra, dec)



class xPeriodicPointSource(xPointSource):

    """Class representing a periodic point source (e.g., a pulsar).

    See :py:class:`ixpeobssim.srcmodel.roi.xCelestialModelComponentBase` for the
    signature of the base class.

    Arguments
    ---------
    ra : float
        The right ascension of the source (in decimal degrees).

    dec : float
        The declination of the source (in decimal degrees).

    ephemeris : :py:class:`ixpeobssim.srcmodel.roi.xEphemeris` object
        The source ephemeris.
    """

    def __init__(self, name, ra, dec, photon_spectrum, polarization_degree,
                 polarization_angle, ephemeris, column_density=0., redshift=0.,
                 identifier=None):
        """Constructor.
        """
        xPointSource.__init__(self, name, ra, dec, photon_spectrum,
                              polarization_degree, polarization_angle,
                              column_density, redshift, identifier)
        self.ephemeris = ephemeris

    def sampling_time_grid(self):
        """Return the phase grid used to sample the lightcurve when generating
        events.
        """
        return numpy.linspace(0., 1., self._count_spectrum_ny)

    def _rvs_seed_event_list(self, parent_roi, irf_set, **kwargs):
        """Overloaded class method for periodic sources.

        There's a fait bit of duplicated code, here, but it's not trivial to do
        much better do to the subtle differences between the roled played by
        time and phase for stationary and periodic sources, respectively.
        """
        # Create the count spectrum.
        phase_grid = self.sampling_time_grid()
        count_spectrum = self.create_count_spectrum(irf_set.aeff, phase_grid, **kwargs)
        # Extract the number of events to be generated. Note we have to multiply
        # the integral of the pulse profile by the duration of the observation.
        start_met, duration = self.parse_time_kwargs(**kwargs)
        num_events = self.poisson(duration * count_spectrum.light_curve.norm())
        logger.info('About to generate %d events...', num_events)
        # Let the ephemeis object throw the phase and time_ arrays.
        phase, time_ = self.ephemeris.rvs(count_spectrum.light_curve, start_met,
            duration, num_events)
        # Apply the GTIs.
        time_, _mask = kwargs.get('gti_list').filter_event_times(time_)
        phase = phase[_mask]
        num_events = len(time_)
        if num_events == 0:
            return xEventList()
        event_list = xEventList(time_, self.identifier)
        # Extract the event energies and perform the pha analysis.
        mc_energy = count_spectrum.rvs(phase)
        mc_pha, mc_pi = irf_set.edisp.pha_analysis(mc_energy)
        # Extract the sky coordinates.
        mc_ra, mc_dec = self.rvs_sky_coordinates(num_events)
        mc_x, mc_y = standard_radec_to_xy(mc_ra, mc_dec, parent_roi.ra, parent_roi.dec)
        # Extract the photoelectron emission directions.
        phi = self._rvs_phi(irf_set.modf, mc_energy, phase, mc_ra, mc_dec)
        # Rotate the phi angles in the gpd reference frame
        detphi = phi_to_detphi(phi, irf_set.du_id, kwargs.get('roll'))
        cols = [mc_energy, mc_pha, mc_pi, mc_ra, mc_dec, mc_x, mc_y, phi, detphi]
        event_list.set_seed_columns(*cols)
        return event_list

    def rvs_photon_list(self, parent_roi, irf_set, **kwargs):
        """Extract a random photon list.
        """
        # Build the custom response at the top of the window.
        aeff_spline, _ = build_tow_response(irf_set)
        # Create the count spectrum.
        phase_grid = self.sampling_time_grid()
        count_spectrum = self.create_count_spectrum(aeff_spline, phase_grid, **kwargs)
        # Extract the number of events to be generated.
        start_met, duration = self.parse_time_kwargs(**kwargs)
        num_events = self.poisson(duration * count_spectrum.light_curve.norm())
        logger.info('About to generate %d photons...', num_events)
        # Let the ephemeris object throw the phase and time_ arrays.
        phase, time_ = self.ephemeris.rvs(count_spectrum.light_curve, start_met,
            duration, num_events)
        # Apply the GTIs.
        time_, _mask = kwargs.get('gti_list').filter_event_times(time_)
        phase = phase[_mask]
        num_events = len(time_)
        if num_events == 0:
            return xPhotonList()
        photon_list = xPhotonList(time_, self.identifier)
        energy = count_spectrum.rvs(phase)
        mc_ra, mc_dec = self.rvs_sky_coordinates(num_events)
        # Smear the sky position with the PSF.
        ra, dec, _, _ = self.convolve_sky_direction(mc_ra, mc_dec, parent_roi, irf_set.psf)
        # Project the sky positions onto the gpd reference frame and apply the
        # dithering effect (if any)
        dither_params = parse_dithering_kwargs(**kwargs)
        roll_angle = kwargs.get('roll')
        ra_pnt, dec_pnt = apply_dithering(time_, parent_roi.ra, parent_roi.dec, dither_params)
        # Project the sky positions onto the gpd reference frame.
        detx, dety = sky_to_gpd(ra, dec, time_, ra_pnt, dec_pnt, irf_set.du_id, roll_angle)
        pol_deg = self.polarization_degree(energy, phase, mc_ra, mc_dec)
        pol_ang = self.polarization_angle(energy, phase, mc_ra, mc_dec)
        # Rotate the polarization angle from the sky reference frame to the
        # GPD reference frame.
        pol_ang = phi_to_detphi(pol_ang, irf_set.du_id, roll_angle)
        photon_list.fill(energy, ra, dec, detx, dety, pol_deg, pol_ang)
        # If the vignetting is enabled, apply it.
        if kwargs.get('vignetting'):
            photon_list.apply_vignetting(irf_set.vign, ra_pnt, dec_pnt)
        return photon_list

    def ephemeris_info(self):
        """Return all the ephemeris information in a form that is suitable for
        display and pretty-printing.
        """
        text = '* Source name: %s\n' % self.name
        text += '* R. A.      : %.6f\n' % self.ra
        text += '* Dec.       : %.6f\n' % self.dec
        met0 = self.ephemeris.met0
        text += '* Epoch MJD  : %.6f (MET %.3f s)\n' % (met_to_mjd(met0), met0)
        text += '* nu0        : %.9f Hz\n' % self.ephemeris.nu0
        text += '* nudot0     : %.6e Hz s^{-1}\n' % self.ephemeris.nudot0
        text += '* nuddot     : %.6e Hz s^{-2}\n' % self.ephemeris.nuddot
        return text

    def __str__(self):
        """String formatting.
        """
        text = xPointSource.__str__(self)
        text += '\n    Ephemeris: %s' % self.ephemeris
        return text



class xBinarySource(xPointSource):

    """Class representing a binary source (e.g., a pulsar in a binary system).

    See :py:class:`ixpeobssim.srcmodel.roi.xModelComponentBase` for the
    signature of the base class.

    Arguments
    ---------
    ra : float
        The right ascension of the source (in decimal degrees).

    dec : float
        The declination of the source (in decimal degrees).

    ephemeris : :py:class:`ixpeobssim.srcmodel.roi.xEphemeris` object
        The source ephemeris.
    """

    def __init__(self, name, ra, dec, photon_spectrum, polarization_degree,
                 polarization_angle, ephemeris, column_density=0., redshift=0.,
                 identifier=None):

        """Constructor.
            """
        xPointSource.__init__(self, name, ra, dec, photon_spectrum, polarization_degree,
                              polarization_angle, column_density, redshift, identifier)
        self.ephemeris = ephemeris

    def _rvs_seed_event_list(self, parent_roi, irf_set, **kwargs):
        """Overloaded class method for periodic sources.

        There's a fait bit of duplicated code, here, but it's not trivial to do
        much better do to the subtle differences between the roled played by
        time and phase for stationary and periodic sources, respectively.
        """
        start_met, duration = self.parse_time_kwargs(**kwargs)
        # Create the underlying count spectrum object.
        phase_grid = numpy.linspace(0., 1., self._count_spectrum_ny)
        count_spectrum = self.create_count_spectrum(irf_set.aeff, phase_grid, **kwargs)
        # Extract the events number to be generated
        average_num_events = duration * count_spectrum.light_curve.norm()
        num_events = numpy.random.poisson(average_num_events)
        logger.info('About to generate %d events...', average_num_events)
        psr_shape = count_spectrum.light_curve
        # Extract times
        time_ = time_list(psr_shape, start_met, self.ephemeris, num_events, duration)
        time_.sort()
        #Extract the event phase
        phase = phase_function(time_, start_met, self.ephemeris.nu(start_met),
            self.ephemeris.nudot(start_met), self.ephemeris.nuddot)
        ph = phase - numpy.floor(phase)
        # Extract the event energies.
        energy = count_spectrum.rvs(time_)
        # Convert into pha and pi
        pha = irf_set.edisp.energy_to_channel(energy)
        pi = pha
        # Extract the sky coordinates
        ra, dec = self.rvs_sky_coordinates(num_events)
        x, y = standard_radec_to_xy(ra, dec, parent_roi.ra, parent_roi.dec)
        # Extract the photoelectron emission directions
        phi = self._rvs_phi(irf_set.modf, energy, ph, ra, dec)
        detphi = phi
        # Build and return the actual event list.
        event_list = xEventList(time_, self.identifier)
        cols = [energy, pha, pi, ra, dec, x, y, phi, detphi]
        event_list.set_seed_columns(*cols)
        return event_list

    def rvs_photon_list(self, parent_roi, irf_set, **kwargs):
        """Placeholder.
        """
        raise NotImplementedError

    def __str__(self):
        """String formatting.
        """
        text = xPointSource.__str__(self)
        text += '\n    Ephemeris: %s' % self.ephemeris
        return text



class xUniformDisk(xCelestialModelComponentBase):

    """Class representing a uniform disk.

    See :py:class:`ixpeobssim.srcmodel.roi.xCelestialModelComponentBase` for the
    signature of the base class.

    Arguments
    ---------
    ra : float
        The right ascension of the disk center (in decimal degrees).

    dec : float
        The declination of the disk center (in decimal degrees).

    radius : float
        The radius of the disk (in degrees).
    """

    def __init__(self, name, ra, dec, radius, photon_spectrum,
                 polarization_degree, polarization_angle, column_density=0.,
                 redshift=0., identifier=None):
        """Constructor.
        """
        args = name, photon_spectrum, polarization_degree, polarization_angle,\
            column_density, redshift, identifier
        xCelestialModelComponentBase.__init__(self, *args)
        self.ra = ra
        self.dec = dec
        self.radius = radius

    # pylint: disable=arguments-differ
    def build_intensity_map(self, wcs_):
        """Overloaded method.
        """
        angsep = angular_separation(self.ra, self.dec, *wcs_to_sky_meshgrid(wcs_))
        mask = angsep <= self.radius
        data = numpy.full(wcs_.array_shape, 0.)
        data[mask] = 1. / mask.sum()
        data = data.transpose()
        return data

    def rvs_sky_coordinates(self, size=1):
        """Generate random coordinates for the model component.

        This is returning an array of the proper length with identical values.

        The algorithm is taken from
        http://mathworld.wolfram.com/DiskPointPicking.html

        Arguments
        ---------
        size : float
            The number of sky coordinate pairs to be generated.
        """
        r = self.radius * numpy.sqrt(numpy.random.sample(size))
        theta = numpy.random.uniform(0, 2. * numpy.pi, size)
        ra = self.ra + (r * numpy.cos(theta) / numpy.cos(numpy.radians(self.dec)))
        dec = self.dec + r * numpy.sin(theta)
        return (ra, dec)

    def __str__(self):
        """String formatting.
        """
        text = xCelestialModelComponentBase.__str__(self)
        text += '\n    Radius: %s deg' % self.radius
        return text



class xGaussianDisk(xCelestialModelComponentBase):

    """Class representing a (azimuthally simmetric) gaussian disk.

    See :py:class:`ixpeobssim.srcmodel.roi.xCelestialModelComponentBase` for the
    signature of the base class.

    Arguments
    ---------
    ra : float
        The right ascension of the disk center (in decimal degrees).

    dec : float
        The declination of the disk center (in decimal degrees).

    sigma : float
        The root mean square of the disk (in degrees).
    """

    def __init__(self, name, ra, dec, sigma, photon_spectrum,
                 polarization_degree, polarization_angle, column_density=0.,
                 redshift=0., identifier=None):
        """Constructor.
        """
        args = name, photon_spectrum, polarization_degree, polarization_angle,\
            column_density, redshift, identifier
        xCelestialModelComponentBase.__init__(self, *args)
        self.ra = ra
        self.dec = dec
        self.sigma = sigma
        self.__mean = [self.ra, self.dec]
        cov00 = (sigma / numpy.cos(numpy.radians(self.dec)))**2.
        cov11 = sigma**2.
        self.__cov = [[cov00, 0.], [0., cov11]]

    # pylint: disable=arguments-differ
    def build_intensity_map(self, wcs_):
        """Overloaded method.
        """
        pdf = scipy.stats.multivariate_normal(self.__mean, self.__cov).pdf
        data = pdf(numpy.dstack(wcs_to_sky_meshgrid(wcs_)))
        data /= data.sum()
        data = data.transpose()
        return data

    def rvs_sky_coordinates(self, size=1):
        """Generate random coordinates for the model component.

        This is returning an array of the proper length with identical values.

        Arguments
        ---------
        size : float
            The number of sky coordinate pairs to be generated.
        """
        rvs = numpy.random.multivariate_normal(self.__mean, self.__cov, size)
        ra, dec = rvs[:, 0], rvs[:, 1]
        return (ra, dec)

    def __str__(self):
        """String formatting.
        """
        text = xCelestialModelComponentBase.__str__(self)
        text += '\n    Sigma: %s deg' % self.sigma
        return text



class xUniformAnnulus(xCelestialModelComponentBase):

    """Class representing a uniform annulus.

    See :py:class:`ixpeobssim.srcmodel.roi.xCelestialModelComponentBase` for the
    signature of the base class.

    Arguments
    ---------
    ra : float
        The right ascension of the annulus center (in decimal degrees).

    dec : float
        The declination of the annulus center (in decimal degrees).

    rmin : float
        The minimium radius of the annulus (in degrees).

    rmax :  float
        The maximium radius of the annulus (in degrees).
    """

    def __init__(self, name, ra, dec, rmin, rmax, photon_spectrum,
                 polarization_degree, polarization_angle, column_density=0.,
                 redshift=0., identifier=None):
        """Constructor.
        """
        args = name, photon_spectrum, polarization_degree, polarization_angle,\
            column_density, redshift, identifier
        xCelestialModelComponentBase.__init__(self, *args)
        self.ra = ra
        self.dec = dec
        assert rmax > rmin
        self.rmin = rmin
        self.rmax = rmax

    # pylint: disable=arguments-differ
    def build_intensity_map(self, wcs_):
        """Overloaded method.
        """
        angsep = angular_separation(self.ra, self.dec, *wcs_to_sky_meshgrid(wcs_))
        mask = numpy.logical_and(angsep >= self.rmin, angsep <= self.rmax)
        data = numpy.full(wcs_.array_shape, 0.)
        data[mask] = 1. / mask.sum()
        data = data.transpose()
        return data

    def rvs_sky_coordinates(self, size=1):
        """Generate random coordinates for the model component.

        This is returning an array of the proper length with identical values.

        The algorithm is taken from
        http://mathworld.wolfram.com/DiskPointPicking.html

        Arguments
        ---------
        size : float
            The number of sky coordinate pairs to be generated.
        """
        r = self.rmin + (self.rmax - self.rmin) * numpy.sqrt(numpy.random.sample(size))
        theta = numpy.random.uniform(0, 2. * numpy.pi, size)
        ra = self.ra + (r * numpy.cos(theta) / numpy.cos(numpy.radians(self.dec)))
        dec = self.dec + r * numpy.sin(theta)
        return (ra, dec)

    def __str__(self):
        """String formatting.
        """
        text = xCelestialModelComponentBase.__str__(self)
        text += '\n    Radius: %s--%s deg' % (self.rmin, self.rmax)
        return text




class xExtendedSource(xCelestialModelComponentBase):

    """Class representing an extended source.

    See :py:class:`ixpeobssim.srcmodel.roi.xCelestialModelComponentBase` for the
    signature of the base class.

    Arguments
    ---------
    img_file_path : string
        The path to the FITS file containing the image of the source.
    """

    def __init__(self, name, img_file_path, photon_spectrum,
                 polarization_degree, polarization_angle, column_density=0.,
                 redshift=0., identifier=None):
        """Constructor.
        """
        args = name, photon_spectrum, polarization_degree, polarization_angle,\
            column_density, redshift, identifier
        xCelestialModelComponentBase.__init__(self, *args)
        self.image = xFITSImage(img_file_path)

    def rvs_sky_coordinates(self, size=1):
        """Generate random coordinates for the model component.

        Arguments
        ---------
        size : float
            The number of sky coordinate pairs to be generated.
        """
        return self.image.rvs_coordinates(size)

    def __str__(self):
        """String formatting.
        """
        text = xCelestialModelComponentBase.__str__(self)
        text += '\n    FITS image: %s' % self.image
        return text



class xChandraObservation(xModelComponentBase):

    """ Class representing a source taken from a Chandra observation.

    Arguments
    ---------
    name : string
        The name of the source.

    polarization_degree : function
        The function object representing the polarization degree.

    polarization_angle : function
        The function object representing the polarization angle.

    region
        The optional region to select the photon list from.

    exclude : bool
        The optional flag to exclude the selected region from the simulation.

    identifier : int
        A unique identifier of the source within a ROI model. Defaults to None
        and is, generally, automatically assigned when building the ROI (i.e.,
        you don't have to worry about it, but it's handy to have in the
        output event list).
    """

    def __init__(self, name, polarization_degree, polarization_angle,
                 region=None, exclude=False, identifier=None):
        """Constructor.
        """
        xModelComponentBase.__init__(self, name, identifier)
        self.polarization_degree = polarization_degree
        self.polarization_angle = polarization_angle
        self.region = region
        self.exclude = exclude

    @staticmethod
    def _time_scaling(scale, energy, ra, dec):
        """Repeat the given arrays according to the scale parameter.

        Warning
        -------
        I have no idea what this is doing. Some refactoring would not hurt?
        """
        num_rep = int(scale[1])
        index = 1 + int(scale[0] * len(energy))
        mc_energy = numpy.append(numpy.tile(energy, num_rep),
                                 energy[:index])
        mc_ra = numpy.append(numpy.tile(ra, num_rep), ra[:index])
        mc_dec = numpy.append(numpy.tile(dec, num_rep), dec[:index])
        return mc_energy, mc_ra, mc_dec

    def _rvs_seed_event_list(self, parent_roi, irf_set, **kwargs):
        """Generate an event list with random time and energy values for the
        model component and for a given effective area and simulation setup.

        Warning
        -------
        The following method works if the effective area ratio is less than 1
        (that is the case of IXPE vs. Chandra).
        """
        start_met, duration = self.parse_time_kwargs(**kwargs)
        roll_angle = kwargs.get('roll')
        mc_energy, mc_ra, mc_dec, mc_effexp =\
            parent_roi.filter_events(self.region, self.exclude)
        num_events = len(mc_energy)
        if num_events == 0:
            return xEventList()
        # Compute for each CHANDRA event the relative exposure ratio and from that
        # derive, using Poisson statistics, the repetion of the event in IXPE
        expect_repeat = irf_set.aeff(mc_energy) * duration / mc_effexp
        # Take care of potential nan or negative values
        expect_repeat[numpy.logical_not(expect_repeat > 0.)] = 0.
        num_repeat = numpy.random.poisson(lam=expect_repeat)
        # Repeat the events
        mc_energy = numpy.repeat(mc_energy, num_repeat)
        mc_ra = numpy.repeat(mc_ra, num_repeat)
        mc_dec = numpy.repeat(mc_dec, num_repeat)
        # Extract times uniformly in the range and initialize the event list.
        time_ = self.uniform_time(len(mc_energy), start_met, duration)
        # Apply the GTIs.
        time_, _mask = kwargs.get('gti_list').filter_event_times(time_)
        # Need to filter the three orignal columns to cope with the events
        # discarded due to the GTI and Earth occultation selections.
        mc_energy = mc_energy[_mask]
        mc_ra = mc_ra[_mask]
        mc_dec = mc_dec[_mask]
        num_events = len(time_)
        if num_events == 0:
            return xEventList()
        event_list = xEventList(time_, self.identifier)
        # Mind all the derived quantities need to be calculated at the very
        # end in order not to have column lenght mismatches.
        mc_pha, mc_pi = irf_set.edisp.pha_analysis(mc_energy)
        mc_x, mc_y = standard_radec_to_xy(mc_ra, mc_dec, parent_roi.ra, parent_roi.dec)
        # Extract the photoelectron emission directions.
        _pd = self.polarization_degree(mc_energy, time_, mc_ra, mc_dec)
        _pa = self.polarization_angle(mc_energy, time_, mc_ra, mc_dec)
        phi = irf_set.modf.rvs_phi(mc_energy, _pd, _pa)
        # Rotate the phi angles in the gpd reference frame
        detphi = phi_to_detphi(phi, irf_set.du_id, roll_angle)
        cols = [mc_energy, mc_pha, mc_pi, mc_ra, mc_dec, mc_x, mc_y, phi, detphi]
        event_list.set_seed_columns(*cols)
        return event_list

    def rvs_event_list(self, parent_roi, irf_set, **kwargs):
        """Extract a random event list for the model component.
        """
        event_list = self._rvs_seed_event_list(parent_roi, irf_set, **kwargs)
        xCelestialModelComponentBase.convolve_event_list(event_list, parent_roi,
                                                         irf_set, **kwargs)
        return event_list

    def rvs_photon_list(self, parent_roi, irf_set, **kwargs):
        """Extract a random photon list.
        """
        start_met, duration = self.parse_time_kwargs(**kwargs)
        # Build the custom response at the top of the window.
        aeff_spline, _ = build_tow_response(irf_set)
        # Extract CHANDRA events
        mc_energy, mc_ra, mc_dec, mc_effexp =\
            parent_roi.filter_events(self.region, self.exclude)
        num_events = len(mc_energy)
        if num_events == 0:
            return xPhotonList()
        # Compute for each CHANDRA event the relative exposure ratio and from that
        # derive, using Poisson statistics, the repetion of the event in IXPE
        expect_repeat = aeff_spline(mc_energy) * duration / mc_effexp
        # Take care of potential nan or negative values
        expect_repeat[numpy.logical_not(expect_repeat > 0.)] = 0.
        num_repeat = numpy.random.poisson(lam=expect_repeat)
        # Repeat the events
        mc_energy = numpy.repeat(mc_energy, num_repeat)
        mc_ra = numpy.repeat(mc_ra, num_repeat)
        mc_dec = numpy.repeat(mc_dec, num_repeat)
        # Extract times uniformly in the range and initialize the event list.
        time_ = self.uniform_time(len(mc_energy), start_met, duration)
        # Apply the GTIs.
        time_, _mask = kwargs.get('gti_list').filter_event_times(time_)
        mc_energy = mc_energy[_mask]
        mc_ra = mc_ra[_mask]
        mc_dec = mc_dec[_mask]
        num_events = len(time_)
        if num_events == 0:
            return xPhotonList()
        photon_list = xPhotonList(time_, self.identifier)
        # Smear the sky position with the PSF.
        ra, dec, _, _ = self.convolve_sky_direction(mc_ra, mc_dec, parent_roi, irf_set.psf)
        # Project the sky positions onto the gpd reference frame and apply the
        # dithering effect (if any)
        dither_params = parse_dithering_kwargs(**kwargs)
        roll_angle = kwargs.get('roll')
        ra_pnt, dec_pnt = apply_dithering(time_, parent_roi.ra, parent_roi.dec, dither_params)
        # Project the sky positions onto the gpd reference frame.
        detx, dety = sky_to_gpd(ra, dec, time_, ra_pnt, dec_pnt, irf_set.du_id, roll_angle)
        pol_deg = self.polarization_degree(mc_energy, time_, mc_ra, mc_dec)
        pol_ang = self.polarization_angle(mc_energy, time_, mc_ra, mc_dec)
        # Rotate the polarization angle from the sky reference frame to the
        # GPD reference frame.
        pol_ang = phi_to_detphi(pol_ang, irf_set.du_id, roll_angle)
        photon_list.fill(mc_energy, ra, dec, detx, dety, pol_deg, pol_ang)
        # If the vignetting is enabled, apply it.
        if kwargs.get('vignetting'):
            photon_list.apply_vignetting(irf_set.vign, ra_pnt, dec_pnt)
        return photon_list



class xROIModel(OrderedDict):

    """Class describing a full ROI (region of interest) model.

    This is essentially an (ordered) collection of component objects
    (i.e., instances of classes inheriting from xCelestialModelComponentBase)
    than can be accessed by source name.

    Arguments
    ---------
    ra_center : float
        The right ascension of the center of the ROI (in decimal degrees).

    dec_center : float
        The declination of the center of the ROI (in decimal degrees).
    """

    def __init__(self, ra_center, dec_center, *sources):
        """Constructor.
        """
        OrderedDict.__init__(self)
        self.ra = ra_center
        self.dec = dec_center
        self.add_sources(*sources)

    def first_source(self):
        """Return the first source (by insertion) in the ROI.
        """
        return self[next(iter(self))]

    # pylint: disable=inconsistent-return-statements
    def source_by_name(self, name=None):
        """Retrieve a source by name.

        If None (or no argument is passed) this returns the first source in the ROI.
        """
        if name is None:
            return self.first_source()
        try:
            return self[name]
        except KeyError:
            abort('ROI model has no source named "%s"' % name)

    # pylint: disable=inconsistent-return-statements
    def source_by_id(self, uid):
        """Retrieve a source by insertion index.
        """
        try:
            return list(self.values())[uid]
        except IndexError:
            abort('ROI model has no source @ id %d' % uid)

    def add_source(self, source):
        """Add a source to the ROI.
        """
        if source.name in self:
            logger.error('ROI already contains source "%s"', source.name)
            abort('Please fix the configuration file')
        source.identifier = len(self)
        self[source.name] = source

    def add_sources(self, *sources):
        """Add an arbitrary number of sources to the ROI.
        """
        for source in sources:
            self.add_source(source)

    def __add__(self, other):
        """Combine different ROI models.
        """
        assert self.__class__.__name__ == other.__class__.__name__
        roi_model = xROIModel(self.ra, self.dec)
        for source in self.values():
            roi_model.add_source(source)
        for source in other.values():
            roi_model.add_source(source)
        return roi_model

    def __str__(self):
        """String formatting.
        """
        txt = 'ROI centered at (%.4f, %.4f):\n' % (self.ra, self.dec)
        for source in self.values():
            txt += '- %s\n' % source
        return txt.strip('\n')

    # pylint: disable=unused-argument, no-self-use
    def _prepare_event_list(self, **kwargs):
        """Initialize an empty event list.

        This is purely a convenience function meant to allow subclasses
        (e.g., xChandraROIModel) to do some initial working before starting
        throwing random numbers.
        """
        return xEventList()

    def _prepare_photon_list(self, **kwargs):
        """Initialize an empty photon list.

        This is purely a convenience function meant to allow subclasses
        (e.g., xChandraROIModel) to do some initial working before starting
        throwing random numbers.

        Note that we really need this segnature for the function, although in
        this base class we are not using any of the arguments. xChandraROIModel
        will. The benefit is that we can factor out all the code in common
        in rvs_event_list(), which subclasses don't have to overload, anymore.
        """
        return xPhotonList()

    def rvs_event_list(self, irf_set, **kwargs):
        """Extract an event list for the full ROI.

        Arguments
        ---------
        irf_set : ixpeobssim.irf.xIRFSet` object.
            The set of instrument response functions to be used.

        Warning
        -------
        The sampling_time should not be the same for all sources, and each
        source should be able to decide its own in a sensible way.
        (See issue #44.)
        """
        event_list = self._prepare_event_list(**kwargs)
        for source in self.values():
            logger.info('Generating event list for "%s"...', source.name)
            event_list += source.rvs_event_list(self, irf_set, **kwargs)
        return event_list

    def rvs_photon_list(self, irf_set, **kwargs):
        """Extract a photon list for the full ROI.

        This was added to support xpphotonlist.py.
        """
        photon_list = self._prepare_photon_list(**kwargs)
        for source in self.values():
            logger.info('Generating photon list for "%s"...', source.name)
            photon_list += source.rvs_photon_list(self, irf_set, **kwargs)
        return photon_list



class xChandraROIModel(xROIModel):

    """Class describing a Chandra ROI (region of interest) model.

    This is essentially an (ordered) collection of component objects
    (i.e., instances of classes inheriting from xCelestialModelComponentBase)
    than can be accessed by source name.

    Arguments
    ---------
    evt_file_path : string
        The path to the FITS file containing the Chandra event list.
    """

    def __init__(self, evt_file_path, acis):
        """Constructor.
        """
        assert acis in ['I', 'S']
        self.acis = acis
        self.evt_file_path = evt_file_path
        logger.info('Reading input Chandra photon list %s...', evt_file_path)
        with fits.open(evt_file_path) as hdu_list:
            hdu_list.info()
            evt_header = hdu_list['EVENTS'].header
            ra_pnt, dec_pnt = chandra.pointing(evt_header)
            xROIModel.__init__(self, ra_pnt, dec_pnt)
            self._wcs(evt_header)
            try:
                self.obs_time = chandra.livetime(evt_header)
            except KeyError:
                self.obs_time = chandra.gti(hdu_list['GTI'].data)
            logger.info('Total Chandra observation time: %f s.', self.obs_time)
            self._load_evt(hdu_list['EVENTS'].data)

    def _wcs(self, evt_header):
        """Retrieve WCS information from the header.
        """
        self.wcs = wcs.WCS(evt_header, keysel=['pixel'], naxis=[5, 6])
        self.wcs.wcs.colax = [0, 0]
        # We need to explicitly set the 'LONPOLE' key because the automatic WCS parsing fails
        # to correctly retrieve it.
        self.wcs.wcs.lonpole = 180.

    def _load_evt(self, evt_data):
        """Read and save the relevant chandra event columns.
        """
        # Convert ev to keV.
        energy_c = 0.001 * evt_data['energy']
        # We cut here in energy to avoid to take the bunch of events with energy
        # greater than 10 keV (probably due to pile-up) and those with energy
        # lower than the minimum energy for IXPE
        _mask = (energy_c > 1.) * (energy_c < 10.)
        self.energy_c = energy_c[_mask]
        x_c = evt_data['x'][_mask]
        y_c = evt_data['y'][_mask]
        self.ra_c, self.dec_c = self.wcs.wcs_pix2world(x_c, y_c, True)
        theta_c = angular_separation(self.ra_c, self.dec_c, self.ra, self.dec)
        self.theta_c = degrees_to_arcmin(theta_c)
        if 'FLUX' in evt_data.columns.names:
            logger.info('FLUX column found in the Chandra event list.')
            logger.info('Using this data to compute the exposure per event...')
            flux_c = erg_to_keV(evt_data['FLUX'][_mask])
            self.effexp_c = self.energy_c / flux_c
        else:
            logger.info('FLUX column NOT found in the Chandra event list.')
            logger.info('Using standard Chandra IRFs to compute the average exposure...')
            self._load_irfs()
            self.effexp_c = self.vign(self.energy_c, self.theta_c) * \
                self.aeff(self.energy_c) * self.obs_time

    def _reset_mask(self, value=False):
        """Set all the elements of filter mask to False (or True).
        """
        assert value in [True, False]
        self.filter_mask = numpy.full(self.energy_c.shape, value, dtype=bool)

    def _check_overlap(self, mask):
        """Return True if the region corresponding to the given mask is
        overlapping one of the others.
        """
        return numpy.logical_and(self.filter_mask, mask).any()

    def _load_irfs(self):
        """Load the Chandra effective area and vignetting.
        """
        detname = 'ACIS-%s' % self.acis
        self.aeff = chandra.load_arf(detname)
        self.vign = chandra.load_vign()

    def filter_events(self, region, exclude=False):
        """Return the filtered event arrays with coordinates inside the given
        region.

        Arguments
        ---------
        region : region instance or None
            The region to select the photon list from (None to take the whole
            remaining area).

        exclude : bool
            Flag to exclude the selected region from the simulation.
        """
        if region is not None:
            _mask = ds9_region_filter_sky(self.ra_c, self.dec_c, self.wcs, region)
            if self._check_overlap(_mask):
                abort('Overlapping region: %s' % region)
            self.filter_mask = numpy.logical_or(_mask, self.filter_mask)
        else:
            _mask = numpy.logical_not(self.filter_mask)
            self._reset_mask(True)
        if exclude:
            return [[] for i in range(4)]
        return [self.energy_c[_mask], self.ra_c[_mask], self.dec_c[_mask],\
                self.effexp_c[_mask]]

    def _prepare_event_list(self, **kwargs):
        """Overloaded method.
        """
        self._reset_mask()
        event_list = xEventList()
        duration = kwargs.get('duration')
        logger.info('Setting the observation time to %d s...', duration)
        return event_list

    def _prepare_photon_list(self, **kwargs):
        """Overloaded method.
        """
        self._reset_mask()
        photon_list = xPhotonList()
        duration = kwargs.get('duration')
        if duration is not None:
            logger.info('Setting the observation time to %d s...', duration)
        return photon_list

    def __str__(self):
        """String formatting.
        """
        text = 'Chandra FITS file: %s' % self.evt_file_path
        text += '\n    %s' % xROIModel.__str__(self)
        return text

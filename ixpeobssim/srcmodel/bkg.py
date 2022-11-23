#!/urs/bin/env python
#
# Copyright (C) 2021--2022, the ixpeobssim team.
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

"""Background model components.
"""

import os

from astropy.io import fits
import numpy

from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.core.rand import xUnivariateGenerator
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.evt.event import xEventList
from ixpeobssim.evt.fmt import standard_radec_to_xy
from ixpeobssim.evt.ixpesim import xPhotonList
from ixpeobssim.irf.ebounds import PI_ENERGY_MIN, PI_ENERGY_MAX, ENERGY_STEP
from ixpeobssim.irfgen.gpd import GPD_FILL_TEMPERATURE, GPD_TYPICAL_ASYMTPTOTIC_PRESSURE
from ixpeobssim.irfgen.gpd import xQeffDataInterface
from ixpeobssim.instrument.gpd import phi_to_detphi
from ixpeobssim.instrument import gpd
from ixpeobssim.instrument.gpd import phi_to_detphi, GPD_PHYSICAL_MAX_RADIUS,\
    GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_Y, within_gpd_physical_area
from ixpeobssim.instrument.mma import gpd_to_sky, parse_dithering_kwargs, apply_dithering
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.srcmodel.roi import xModelComponentBase
from ixpeobssim.srcmodel.roi import xUniformDisk
from ixpeobssim.srcmodel.spectrum import power_law, cutoff_power_law, \
    xSourceSpectrum, load_spectral_spline
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.units_ import arcmin_to_degrees


# pylint: disable=invalid-name, too-many-locals



class xRadialBackgroundGenerator(xUnivariateGenerator):

    """Univariate generator sampling the radial coordinate on the detector surface
    for a background component whose radial distribution, normalized by the area,
    depends linearly on the radius.

    This was introduced in https://github.com/lucabaldini/ixpeobssim/issues/663

    The need for an instrumental background component that is not uniform in
    detector coordinates stems from the analysis of a number of observations.
    Interestingly enough, the radial slope of the background counts per unit
    area seems to be different for different observation.

    .. warning::

       While initially I was hoping to code this by sampling two independent
       random variables, within the proper bounds, on the x and y coordinates,
       it turned out that an azimuthally symmetric bivariate pdf cannot be
       expressed as the product of two independent variables on a square, and
       we had to resort to sampling r over a circle and trimming in the fiducial
       rectangle after the fact. This is hugly, as we need to guess in advance
       how much random numbers we have to throw so that we end up with enough
       counts after the trimming, but so life goes.

    This is essentially an univariate generator whose pdf is a slight generalization
    of the function f(r) = 2r (for r = 1) that one would use to throw random
    numbers uniformly distributed in a circle:

    .. math::

       p(r) = \\left( 1 - \\frac{\\alpha}{2}\\right) \\frac{r}{h} +
       \\alpha \\left(\\frac{r}{h}\\right)^2

    The radial slope alpha repesents the fractional half-excursion of the variation
    across the size h of the fiducial rectangle. For alpha = 0 the detector
    position are distributed uniformly over the fiducial rectangle. For alpha = 2
    the radial dependence is maximal, and the density of events is zero at the
    center of the detector. For alpha = -1 (and assuming a square fiducial region)
    the pdf approaches zero at the boundary of the circle.

    Arguments
    ---------
    radial_slope : float
        The slope of the radial profile, that is, the fractional half-excursion
        of the variation across the size of the fiducial rectangle.
    """

    HALF_SIDE = 0.5 * (GPD_PHYSICAL_HALF_SIDE_X + GPD_PHYSICAL_HALF_SIDE_Y)

    def __init__(self, radial_slope, num_points=100):
        """Constructor.
        """
        if radial_slope > 2. or radial_slope < -1.:
            raise RuntimeError('Invalid background radial slope (%.3f)' % radial_slope)
        self.radial_slope = radial_slope
        r = numpy.linspace(0., GPD_PHYSICAL_MAX_RADIUS, num_points)
        xUnivariateGenerator.__init__(self, r, self.pdf(r, self.radial_slope))

    @staticmethod
    def pdf(r, radial_slope):
        """Small function encapsulating the underlying pdf for the random generator.

        Note we are using the average physical size of the GPD along the x and
        y directions as an effective value for the radial parametrization.
        """
        r = r / xRadialBackgroundGenerator.HALF_SIDE
        return (1. - radial_slope / 2.) * r  + radial_slope * r**2.

    @staticmethod
    def polar_to_cartesian(r, phi):
        """Convert an array of polar coordinates in the plane into the corresponding
        cartesian coordinates.
        """
        x = r * numpy.cos(phi)
        y = r * numpy.sin(phi)
        return x, y

    @staticmethod
    def average_oversample_fraction(radial_slope):
        """Return an heuristic for the average oversample fraction for a given
        radial slope of the profile.

        This is a purely geometric factor that is easy to calculate for a flat
        profile (slope = 0), in which case it reads pi/2 but not trivial in the
        general case. Our approach is to generate events within the smallest circle
        ciscumscribed to the fiducial rectangle on a grid of radial slope values
        and measure the fraction that ends up in the fiducial rectangle itself.
        For completeness, this is calculated in tests/test_radial_bkg.py.

        Note that this is accurate within a few % in the limit of infinite statistics.
        """
        return 1.56258183 + 0.29704984 * radial_slope - 0.05847132 * radial_slope**2.

    def oversampled_size(self, size, radial_slope, safety_factor=1.2, min_size=1000):
        """Return the oversampled size for a given target size and radial slope.

        This is essentially the function that determines how many random numbers
        we need to throw to be sure that, after trimming to the fiducial region,
        we end up to enough events. This is purely heuristic and is based on the
        average_oversample_fraction() function when the statistics is large enough,
        with a minimum bound to be sure we are not killed by statistical fluctuation
        in the small number regime.
        """
        size = safety_factor * size * self.average_oversample_fraction(radial_slope)
        size = round(size)
        return max(size, min_size)

    def rvs_xy(self, size):
        """Extract a set of arrays of coordinate detectors.
        """
        oversize = self.oversampled_size(size, self.radial_slope)
        logger.info('Sampling the background radial profile...')
        logger.info('Initial counts: %d, target counts %d', oversize, size)
        r = self.rvs(oversize)
        phi = 2. * numpy.pi * numpy.random.random(oversize)
        x, y = self.polar_to_cartesian(r, phi)
        mask = within_gpd_physical_area(x, y)
        x = x[mask]
        y = y[mask]
        if len(x) < size:
            raise RuntimeError('Not enough background counts after trimming.')
        x = x[:size]
        y = y[:size]
        return x, y



class xInstrumentalBkg(xModelComponentBase):
    """
    Base class for the instrumental background

    In addition to the base class xModelComponentBase, this class has the
    spectrum member to fully describe the spectrum of the instrumental
    background.

    This essentially corresponds to a uniform, unpolarized source in the GPD
    reference frame, whose events are not convolved with any of the instrument
    response functions, except for the energy dispersion, if explicitely requested.

    Arguments
    ---------
    name : str
        The name of the component

    photon_spectrum : callable
        The photon spectrum in photons cm-2/s-1/keV-1 for the bkg source.

    radial_slope : float
        The radial slope for the distribution in detector coordinates.

    identifier : int, optional
        The source identifier.

    convolve_energy : bool
        Convolve the source spectrum with the energy dispersion (default is False).
    """

    def __init__(self, name, photon_spectrum, radial_slope=0., identifier=None,
                 convolve_energy=False):
        """Constructor.
        """
        xModelComponentBase.__init__(self, name, identifier)
        self.photon_spectrum = photon_spectrum
        self.radial_slope = radial_slope
        self._convolve_energy = convolve_energy

    def _energy_grid(self, irf_set, num_points=250, **kwargs):
        """Return the proper energy grid for the energy convolution.

        This requires some care, as the implications are slightly different
        depending we are convolving the input spectrum with the energy dispersion
        or not.

        If we are convolving the energy, we are doing exactly the same thing
        that we are doing for the celestial sources. If that's not the case
        we use the full 0.04--15 keV range. (Note in this case we're not starting
        from 0 keV not to run into problems with power laws diverging in the
        origin.)
        """
        if self._convolve_energy:
            emin = kwargs.get('emin', irf_set.aeff.xmin())
            emax = kwargs.get('emax', irf_set.aeff.xmax())
        else:
            emin = PI_ENERGY_MIN + ENERGY_STEP
            emax = PI_ENERGY_MAX
        return numpy.linspace(emin, emax, num_points)

    def _seed_columns(self, irf_set, gpd_qeff_correction=False, **kwargs):
        """Fill the seed columns.
        """
        start_met, duration = self.parse_time_kwargs(**kwargs)
        energy_grid = self._energy_grid(irf_set)
        # Note: this is hard-coded and should be changed.
        time_grid = numpy.linspace(start_met, start_met + duration, 100)
        # Create the source spectrum
        if gpd_qeff_correction:
            # If we are using this for a photon list, then we have to de-correct
            # for the effect of the GPD quantum efficiency, as ixpesim, at this
            # stage, only simulates photons at the top of the window. Note this
            # will never be exact, as the instrumental background is a mixture of
            # different things. We might want to simulate the correct misture of
            # charged particles, at some point, but for the time being this will
            # get at least in the right ball-park as far as the number of events
            # is concerned.
            #
            # First thing first, we have to override the default energy grid in
            # order to avoid generating an enormous number of events below 1 keV,
            # where the quantum efficiency is small. Note that, in this case,
            # the bounds are essentially put by hands. It should also be
            # emphasize that, in the ixpesim workflow, there is essentially no
            # way to disengage the energy dispersion, and we will probably never
            # get things completely right.
            logger.info('Overriding the default energy grid for the efficiency decorrection...')
            energy_grid = numpy.linspace(kwargs.get('emin', 1.), PI_ENERGY_MAX, 250)
            # Then we have to retrieve the proper quantum efficiency data
            # for the decorrection to be applied.
            qeff_data = xQeffDataInterface()
            # And for this we need the right temperature and asymptotic pressure.
            # It is unfortunate that this information is buried into the header
            # comments of the response files. As in ixpeobssim.evt.ixpesim, we
            # try and start from sensible defaults...
            temperature = GPD_FILL_TEMPERATURE
            pressure = GPD_TYPICAL_ASYMTPTOTIC_PRESSURE
            # ... and then we try and do some magic with the header comments.
            for line in irf_set.aeff.header_comments():
                if line.startswith('GPD filling temperature'):
                    temperature = float(line.split(':')[-1].replace('deg C', '').strip())
                    logger.info('GPD filling temperature read from comments: %s', temperature)
                if line.startswith('GPD DME asymptotic pressure'):
                    pressure = float(line.split(':')[-1].replace('mbar', '').strip())
                    logger.info('GPD pressure read from comments: %s', pressure)
            logger.info('Decorrecting the GPD efficiency at %.1f deg C, %.1f mbar',
                temperature, pressure)

            def spectrum(E, t=None):
                """Small nested function to apply the quantum efficiency decorrection
                for the simulation of photon lists.

                Note that, by the time this is called, E and t are both two-dimensional
                arrays on the proper grid for the bivariate spline underlying the
                source spectrum to be build. Since the quantum efficiency is
                purely one dimensional, we do have to take one slice in energy,
                calculate the efficiencty, tile the output on the time grids and
                transpose the thing to have the correct output shape.
                """
                energy = E[:,0]
                qeff = qeff_data.quantum_efficiency(energy, temperature, pressure, contaminants=None)
                qeff = numpy.tile(qeff, (len(time_grid), 1)).T
                return self.photon_spectrum(E, t) / qeff

        else:
            # For a standard event list we do the simple thing.
            spectrum = self.photon_spectrum
        source_spectrum = xSourceSpectrum(energy_grid, time_grid, spectrum)
        # Extract the number of events to be generated.
        num_events = self.poisson(source_spectrum.build_light_curve().norm())
        # Multiply by the total number of events by the detector area.
        # Note the fiducial area must be converted from mm2 to cm2
        num_events = int(num_events * gpd.GPD_PHYSICAL_AREA / 100. + 0.5)
        logger.info('About to generate %d events...', num_events)
        # Extract the event time, sort the values and initialize the event list.
        time_ = source_spectrum.build_light_curve().rvs(num_events)
        time_.sort()
        # Apply the GTIs.
        time_, _ = kwargs.get('gti_list').filter_event_times(time_)
        # Extract the event energies.
        mc_energy = source_spectrum.rvs(time_)
        # Extract the event positions---note this was changed in response to
        # https://github.com/lucabaldini/ixpeobssim/issues/663
        # And we still need to handle the fiducial rectangle properly.
        detx, dety = xRadialBackgroundGenerator(self.radial_slope).rvs_xy(len(time_))
        return time_, mc_energy, detx, dety

    def rvs_event_list(self, parent_roi, irf_set, **kwargs):
        """Overloaded method.
        """
        time_, mc_energy, detx, dety = self._seed_columns(irf_set, **kwargs)
        num_events = len(time_)
        if num_events == 0:
            return xEventList()
        event_list = xEventList(time_, self.identifier)
        event_list.set_detector_position_columns(detx, dety)
        # Perform the pha analysis.
        mc_pha, mc_pi = irf_set.edisp.pha_analysis(mc_energy)
        event_list.set_mc_energy_columns(mc_energy, mc_pha, mc_pi)
        if self._convolve_energy:
            energy, pha, pi = irf_set.edisp.convolve_energy(mc_energy)
        else:
            energy, pha, pi = mc_energy, mc_pha, mc_pi
        event_list.set_energy_columns(energy, pha, pi)
        # ... and azimuthal angle.
        roll_angle = kwargs.get('roll')
        phi = self.uniform_phi(num_events)
        detphi = phi_to_detphi(phi, irf_set.du_id, roll_angle)
        event_list.set_phi_columns(phi, detphi)
        # Apply the dithering effect to the pointing direction (if needed).
        dither_params = parse_dithering_kwargs(**kwargs)
        ra_pnt, dec_pnt = apply_dithering(time_, parent_roi.ra, parent_roi.dec, dither_params)
        ra, dec = gpd_to_sky(detx, dety, time_, ra_pnt, dec_pnt, irf_set.du_id, roll_angle)
        # Sky positions, no difference between measured and mc because
        # we are not convolving with the psf
        x, y = standard_radec_to_xy(ra, dec, parent_roi.ra, parent_roi.dec)
        cols = [mc_energy, mc_pha, mc_pi, ra, dec, x, y, phi, detphi]
        event_list.set_seed_columns(*cols)
        cols = [pha, pi, energy, ra, dec, x, y, detx, dety]
        event_list.set_rec_columns(*cols)
        return event_list

    def rvs_photon_list(self, parent_roi, irf_set, **kwargs):
        """Extract a random photon list.
        """
        time_, mc_energy, detx, dety = self._seed_columns(irf_set, True, **kwargs)
        num_events = len(time_)
        if num_events == 0:
            return xPhotonList()
        photon_list = xPhotonList(time_, self.identifier)
        # Apply the dithering effect to the pointing direction (if needed).
        dither_params = parse_dithering_kwargs(**kwargs)
        roll_angle = kwargs.get('roll')
        ra_pnt, dec_pnt = apply_dithering(time_, parent_roi.ra, parent_roi.dec, dither_params)
        # The instrumental background is unpolarized.
        # Shall we allow for something to be set in the constructor?
        pol_deg = numpy.full(num_events, 0.)
        pol_ang = numpy.full(num_events, 0.)
        pol_ang = phi_to_detphi(pol_ang, irf_set.du_id, roll_angle)
        photon_list.fill(mc_energy, ra_pnt, dec_pnt, detx, dety, pol_deg, pol_ang)
        return photon_list



class xPowerLawInstrumentalBkg(xInstrumentalBkg):

    """Specialized instrumental background source with a power-law spectrum.

    The default values for the spectral parameters are set after Bunner et al.
    1978 (1978ApJ...220..261B), where the authors provide the non X-ray
    background rates for their three detectors. We are using values for the Neon
    filled detector in Table 3 of the paper. We fit the three points
    (energy bins: 0.76--1.6, 1.6--3.0, 3.0--6.0) with a power law. The best-fit
    values for the index and normalization of this power-law are 1.0 and 4.e-4,
    respectively.
    """

    def __init__(self, norm=4.e-4, index=1.0, radial_slope=0.):
        """Constructor.
        """
        xInstrumentalBkg.__init__(self, 'Instrumental background', power_law(norm, index),
            radial_slope)
        self.norm = norm
        self.index = index



class xTemplateInstrumentalBkg(xInstrumentalBkg):

    """Instrumental background based on a template.
    """

    DEFAULT_PATH = os.path.join(IXPEOBSSIM_SRCMODEL, 'ascii', 'bkg_smcx1_01903701.txt')

    def __init__(self, file_path=DEFAULT_PATH, emin=0.1, emax=15., k=1, radial_slope=0.):
        """Constructor.
        """
        self.spline = load_spectral_spline(file_path, emin, emax, k=k)
        spec = lambda E, t=None: self.spline(E)
        xInstrumentalBkg.__init__(self, 'Instrumental background', spec, radial_slope)



class xCelestialBkgBase(xUniformDisk):

    """Base class for the celestial background components.

    Celestial background is assumed to be uniform within the field of view, and
    is therefore modeled as a uniform disk with a radius larger than the
    field of view itself---the normalization should be computed self-consistently
    for the full disk, and when the GPD fiducial cut is applied the photons that
    do not hit the detector are simply thrown away.

    (Other than the fact that the normalization is normalized to the solid angle
    of the component, this is essentially identical to any other celestial
    component.)

    Note that, by construction, the component is unpolarized, and the column
    density and redshift are identically zero.
    """

    _RADIUS_ARCMIN = 9.
    _RADIUS_DEG = arcmin_to_degrees(_RADIUS_ARCMIN)

    def __init__(self, name, ra, dec, photon_spectrum, identifier=None):
        """Constructor.
        """
        kwargs = dict(polarization_degree=constant(0.), polarization_angle=constant(0.),
                      column_density=0., redshift=0., identifier=identifier)
        xUniformDisk.__init__(self, name, ra, dec, self._RADIUS_DEG, photon_spectrum, **kwargs)

    @staticmethod
    def solid_angle(radius):
        """Return the solid angle subtended by a cone of aperture radius.

        See https://en.wikipedia.org/wiki/Solid_angle

        Args
        ----
        radius : float
            The half-apex of the cone in decimal degrees.
        """
        return 2. * numpy.pi * (1. - numpy.cos(numpy.radians(radius)))



class xExtragalacticBkg(xCelestialBkgBase):

    """Derived class for the extragalactic X-ray background.

    All the basic formalism is taken from Gruber et al., 1999
    https://iopscience.iop.org/article/10.1086/307450/pdf

    Essentially in our energy band we model the extragalactic background as
    a broken power law with index 1.29 and break energy 41.13 keV.

    Note 1.29 is the index of the photon spectrum, while the 0.29 in the paper
    refers to the energy spectrum, see
    https://bitbucket.org/ixpesw/ixpeobssim/issues/544
    """

    _SPEC_NORM = 7.877
    _SPEC_INDEX = 1.29
    _SPEC_CUTOFF = 41.13

    def __init__(self, ra, dec, identifier=None):
        """Constructor.
        """
        norm = self._SPEC_NORM * self.solid_angle(self._RADIUS_DEG)
        photon_spectrum = cutoff_power_law(norm, self._SPEC_INDEX, self._SPEC_CUTOFF)
        args = ra, dec, photon_spectrum, identifier
        xCelestialBkgBase.__init__(self, 'Extragalactic background', *args)



class xRosatPSPCResponseMatrix:

    """Simple interface to the ROSAT PSPC response matrix.

    The FITS files was downloaded from the ROSAT HEASARC page, and the
    content seems qualitatively in agreement with the expectations based on the
    XRT collecting area
    https://heasarc.gsfc.nasa.gov/docs/rosat/ruh/handbook/node39.html#figXMAareas
    and the PSPC quantum efficiency
    https://heasarc.gsfc.nasa.gov/docs/rosat/ruh/handbook/node56.html#SECTION00730000000000000000

    The on-axis effective area is included in the response matrix, and to
    extract it we just sum the (un-normalized) pdf in each energy bin.

    The effective area is stored internally as an interpolated spline of order k=1.
    """

    _FILE_PATH = os.path.join(IXPEOBSSIM_SRCMODEL, 'fits', 'pspcc_gain1_256.fits')
    # Energy bounds (in keV) for the R7 channel.
    R7_EMIN = 1.05
    R7_EMAX = 2.04

    def __init__(self):
        """Constructor.
        """
        with fits.open(self._FILE_PATH) as hdu_list:
            data = hdu_list['SPECRESP MATRIX'].data
            energy = 0.5 * (data['ENERG_LO'] + data['ENERG_HI'])
            aeff = data['MATRIX'].sum(axis=1)
        fmt = dict(xlabel='Energy [keV]', ylabel='On-axis effective area [cm$^2$]')
        self.aeff = xInterpolatedUnivariateSpline(energy, aeff, k=1, **fmt)



class xGalacticBkg(xCelestialBkgBase):

    """Derived class for the galactic X-ray background.

    A good reference for this is Tanaka, 2002:
    https://www.aanda.org/articles/aa/pdf/2002/06/aah3204.pdf

    Operationally, I think the most convenient mean to gauge the intensity of the
    Galactic background is the HEASARC background tool
    https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xraybg/xraybg.pl
    and, particularly, the ROSAT count rate in the 1--2 keV energy band (R7).
    """

    _SPEC_INDEX = 1.3
    _NORM_SCALE = 5.e-6

    def __init__(self, ra, dec, rosat_r7_bkg_rate, identifier=None):
        """Constructor.

        Args
        ----
        ra : float
            RA for the target position in decimal degrees

        dec : float
            DEC for the targer position in decimal degrees

        rosat_r7_bkg_rate : float
            ROSAT X-ray background average count rate in the R7 band [1e-6 counts/sec/arcmin^2]
        """
        norm = self._NORM_SCALE * rosat_r7_bkg_rate
        photon_spectrum = power_law(norm, self._SPEC_INDEX)
        args = ra, dec, photon_spectrum, identifier
        xCelestialBkgBase.__init__(self, 'Galactic background', *args)

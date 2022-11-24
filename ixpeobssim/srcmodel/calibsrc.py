#!/urs/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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

"""In-flight calibration sources and calibration-related facilities.
"""

from __future__ import print_function, division

import os

import numpy
from astropy.io import fits

from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.core.hist import xHistogram2d
from ixpeobssim.evt.event import xEventList
from ixpeobssim.instrument.gpd import gpd_map_binning, GPD_PHYSICAL_HALF_SIDE_X,\
    GPD_PHYSICAL_HALF_SIDE_Y
from ixpeobssim.srcmodel.roi import xModelComponentBase
from ixpeobssim.utils.logging_ import logger


# pylint: disable=invalid-name

class xCalibrationSourceBase(xModelComponentBase):

    """Base class for all (on-ground and on-orbit) calibration sources.

    The common denominator to all the calibration sources is:

    * the event rate is constant in time and the counting statistics is Poisson;
    * a specific identifier (>=100) is assigned to each source, and set once and
      for all, independently on the ROI machinery (this is different from the
      celestial sources, where the identifier is typically assigned depending on
      the order in which the ROI is populated);
    * there is no concept of parent ROI, as the calibration sources are supposed
      to be operated one at a time, and independently from the satellite pointing;
    * there is no concept of GTI, the start_met and duration being literally
      all we care about when it comes to decide whether the source is
      activated or not.
    """

    NAME = 'Generic calibration source'
    IDENTIFIER = 100

    def __init__(self, rate):
        """Constructor.
        """
        super().__init__(self.NAME, self.IDENTIFIER)
        self.rate = rate

    def rvs_time(self, start_met, duration):
        """Extract the event times.

        Here we essentially throw a random number with a Poisson distribution
        based on the expected mean of events, and then extract the event times
        at random---uniformly distributed between the start and stop met.

        This is supposed to be the common paradigm for all the calibration sources.
        """
        num_events = self.poisson(duration * self.rate)
        return self.uniform_time(num_events, start_met, duration)

    def rvs_energy(self, size):
        """Extract the event energies.

        Do-nothing method to be re-implemented in derived classes.
        """
        raise NotImplementedError

    def rvs_detxy(self, size):
        """Extract the event positions in detector coordinates.

        Do-nothing method to be re-implemented in derived classes.
        """
        raise NotImplementedError

    def rvs_detphi(self, size):
        """Extract the photoelectron emission direction in detector coordinates.

        Do-nothing method to be re-implemented in derived classes.
        """
        raise NotImplementedError

    # pylint: disable=unused-argument
    def rvs_event_list(self, irf_set, **kwargs):
        """Main interface to generate random events from the calibration source.

        Note that the signature here is different from the celestial sources, as
        we don'y need the reference to the parent ROI.

        Also, since the calibration sources are not meant to be combined, the
        deadtime and the trigged id are calculated right away in the body function.
        """
        start_met = kwargs.get('start_met')
        duration = kwargs.get('duration')
        deadtime = kwargs.get('deadtime')
        stop_met = start_met + duration
        # Extract the event time.
        time_ = self.rvs_time(start_met, duration)
        num_events = len(time_)
        logger.info('About to generate %d event(s)...', num_events)
        if num_events == 0:
            return xEventList()
        event_list = xEventList(time_, self.identifier)
        # Sky positions (Monte Carlo and measured). Note that this has no real
        # meaning for calibration sources, and we do set everything to zero.
        ra = numpy.full(num_events, 0.)
        dec = numpy.full(num_events, 0.)
        x = numpy.full(num_events, 0.)
        y = numpy.full(num_events, 0.)
        event_list.set_mc_sky_position_columns(ra, dec, x, y)
        event_list.set_sky_position_columns(ra, dec, x, y)
        # Detector positions.
        detx, dety = self.rvs_detxy(num_events)
        event_list.set_detector_position_columns(detx, dety)
        # Event energy.
        mc_energy = self.rvs_energy(num_events)
        mc_pha, mc_pi = irf_set.edisp.pha_analysis(mc_energy)
        event_list.set_mc_energy_columns(mc_energy, mc_pha, mc_pi)
        energy, pha, pi = irf_set.edisp.convolve_energy(mc_energy)
        event_list.set_energy_columns(energy, pha, pi)
        # Azimuthal angle. Note that, for what it's worth, the two are the
        # same in detector and sky coordinates.
        detphi = self.rvs_detphi(num_events)
        phi = detphi
        event_list.set_phi_columns(phi, detphi)
        return event_list

    def __str__(self):
        """String formatting.
        """
        return 'Calibration source "%s" (id = %d) @ %.3f Hz average rate' %\
            (self.name, self.identifier, self.rate)



class xMonochromaticUnpolarizedCalibrationSourceBase(xCalibrationSourceBase):

    """Specialized base class for monochromatic unpolarized calibration sources.
    """

    FE_55_KA1_ENERGY = 5.89875

    def __init__(self, rate, energy):
        """Constructor.
        """
        super().__init__(rate)
        self.energy = energy

    def rvs_energy(self, size):
        """Overloade method.
        """
        return numpy.full(size, self.energy)

    def rvs_detxy(self, size):
        """Do-nothing method to be re-implemented in derived classes.
        """
        raise NotImplementedError

    def rvs_detphi(self, size):
        """Overloade method.
        """
        return self.uniform_phi(size)



class xCalibrationSourceImage:

    """Small convenience class descring the morphology of a calibration source
    in detector coordinates.

    .. warning::

        This functionality has significant overlap with the xFITSImage class in
        srcmodel.roi and the two should be refactored. Also: the very concept of
        a map in detector coordinates stored in the form of a FITS file
        has a lot in common with the stuff on the spurious modulation branch, and
        all of these things should be accommodated in a consistent fashion.

    The HALF_SIDE top-level class member was added in response to issue
    https://github.com/lucabaldini/ixpeobssim/issues/668
    and was set to the old value of the fiducial half-side. This is now
    self-contained to the calibration source machinery, and we might have to
    generalize the code to arbitrary fiducial cuts in the future, if we ever want
    to heavily use this functionality (which has not been the case, up to now).
    """

    HALF_SIDE = 7.350

    def __init__(self, file_path):
        """Constructor.
        """
        with fits.open(file_path) as hdu_list:
            data = hdu_list['PRIMARY'].data
            self.data = data
            self.shape = self.data.shape
            # Note the cast to float is a terrible workaround for issue
            # https://bitbucket.org/ixpesw/ixpeobssim/issues/608
            # triggered by numpy 1.22.0
            self.cdf = numpy.cumsum(data.ravel().astype(float))
            self.cdf /= self.cdf[-1]
        self.xbinning, self.ybinning = gpd_map_binning(self.HALF_SIDE, self.HALF_SIDE,
            *self.shape)

    def plot(self):
        """Plot the source image.
        """
        h = xHistogram2d(self.xbinning, self.ybinning, xlabel='x [mm]', ylabel='y [mm]]')
        h.set_content(self.data)
        h.plot()

    def rvs_coordinates(self, size=1):
        """Generate random coordinates based on the image map.
        """
        u = numpy.random.rand(size)
        pixel = numpy.searchsorted(self.cdf, u)
        row, col = numpy.unravel_index(pixel, self.shape)
        dx = self.xbinning[1] - self.xbinning[0]
        dy = self.ybinning[1] - self.ybinning[0]
        x = self.xbinning[row] + numpy.random.uniform(0., dx, size)
        y = self.ybinning[col] + numpy.random.uniform(0., dy, size)
        return x, y



class xCalC(xMonochromaticUnpolarizedCalibrationSourceBase):

    """Class describing the "Cal C" onboard calibration source.
    """

    NAME = 'Cal C'
    IDENTIFIER = 101

    def __init__(self, rate):
        """Constructor.
        """
        super().__init__(rate, self.FE_55_KA1_ENERGY)
        file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'fits', 'onboard_cal_C_img.fits')
        self.image = xCalibrationSourceImage(file_path)

    def rvs_detxy(self, size):
        """Overloaded method.
        """
        return self.image.rvs_coordinates(size)



class xMonochromaticUnpolarizedFlatField(xMonochromaticUnpolarizedCalibrationSourceBase):

    """Unpolarized square flat field at a given energy.
    """

    NAME = 'Flat field'

    def rvs_detxy(self, size):
        """Overloaded method.
        """
        return self.uniform_rectangle(size, GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_Y)



class xCalibrationROIModel(dict):

    """Class describing a custom ROI (region of interest) model for calibration
    sources.

    In contrast to celestial ROI models, calibration ROI models only handle
    one source.
    """

    def __init__(self, source):
        """Constructor.
        """
        self.source = source
        dict.__init__(self, {source.identifier: self.source})
        self.ra = 'N/A'
        self.dec = 'N/A'

    def __str__(self):
        """String formatting.
        """
        return 'Calibration ROI model\n%s' % self.source

    def rvs_event_list(self, irf_set, **kwargs):
        """Extract an event list for the full ROI.
        """
        return self.source.rvs_event_list(irf_set, **kwargs)

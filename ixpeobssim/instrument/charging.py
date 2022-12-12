#!/urs/bin/env python
#
# Copyright (C) 2019, the ixpeobssim team.
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

"""Charging-related facilities.
"""

from __future__ import print_function, division

import numpy
from astropy.io import fits

from ixpeobssim.core.fitsio import xPrimaryHDU, xBinTableHDUBase
from ixpeobssim.core.hist import xGpdMap3d, xHistogramBase, xHistogram2d
from ixpeobssim.instrument.gpd import GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_Y
from ixpeobssim.utils.logging_ import logger, abort


# pylint: disable=invalid-name, no-member, too-many-arguments, too-many-ancestors
# pylint: disable=too-many-locals


DEFAULT_K_C = 1.2e5
DEFAULT_TAU_D = 1.1e5
DEFAULT_DELTA_MAX = 0.067


def charging_tau(energy_flux, k_c=DEFAULT_K_C):
    """Return the time constant for the GEM charging in absence of discharge.
    """
    return k_c / energy_flux


def effective_tau(energy_flux, k_c=DEFAULT_K_C, tau_d=DEFAULT_TAU_D):
    """Return the effective time constant for the GEM in an arbitrary state.
    """
    return 1. / (1. / charging_tau(energy_flux, k_c) + 1. / tau_d)


def asymptotic_gain_delta(energy_flux, k_c=DEFAULT_K_C, tau_d=DEFAULT_TAU_D,
                          delta_max=DEFAULT_DELTA_MAX):
    """Return the asymptotic gain decrease for the GEM in an arbitrary state.
    """
    return -delta_max / (1. + charging_tau(energy_flux, k_c) / tau_d)


def gain_profile(energy_flux, k_c=DEFAULT_K_C, tau_d=DEFAULT_TAU_D):
    """Return a lambda to generate a gain profile, given some model parameters.

    to plot the actual time-profile of the gain.
    """
    tau = effective_tau(energy_flux, k_c, tau_d)
    delta = asymptotic_gain_delta(energy_flux, k_c, tau_d)
    return lambda t: 1. + delta * (1. - numpy.exp(-t / tau))


def calculate_gain(time_, energy_flux, initial_gain=1., k_c=DEFAULT_K_C,
                   tau_d=DEFAULT_TAU_D, delta_max=DEFAULT_DELTA_MAX):
    """Calculate the gain as a function of time for a given input energy flux.

    The time grid and the corresponding energy flux are passed as numpy arrays
    (of the same shape). The system is integrated using finite differences.
    """
    gain = [initial_gain]
    for dt, r in zip(numpy.diff(time_), energy_flux[1:]):
        g = gain[-1]
        delta = 1. - g
        g -= r * delta_max / k_c * (1 - delta / delta_max) * dt
        g += delta / tau_d * dt
        gain.append(g)
    gain = numpy.array(gain)
    return gain



class xChargingPrimaryHDU(xPrimaryHDU):

    """Charging map file primary header.
    """

    HEADER_KEYWORDS = [
        ('TELESCOP', 'IXPE'     , 'Telescope name'),
        ('INSTRUME', 'GPD'      , 'Instrument name'),
        ('DETNAM'  , 'N/A'      , 'Name of the logical detector unit'),
        ('DET_ID'  , 'N/A'      , 'Name of the physical detector unit'),
        ('CREATOR' , 'N/A'      , 'Creator app'),
        ('ORIGIN'  , 'IXPE team', 'Organization responsible for the data'),
        ('DATE'    , 'N/A'      , 'File creation date')
    ]
    HEADER_COMMENTS = []



class xBinTableHDUCharging(xBinTableHDUBase):

    """Binary table description for the CHRG_MAP extension of the observation
    output files.
    """

    NAME = 'CHRG_MAP'
    HEADER_KEYWORDS = [
        ('VERSION' , 1    ,  'Extension version number'),
        ('CVSD0001', 'N/A',  'Date when this file should first be used'),
        ('CVST0001', 'N/A',  'Time of day when this file should first be used'),
        ('NUM_BINS', 'N/A',  'Number of bins per side of the map')
    ]
    HEADER_COMMENTS = ['This extension provides a map of the detector '\
                'charging status expressed as a fraction of its maximum value.']
    DATA_SPECS = [
        ('BINX', 'I', None, 'Index for the x coordinate'),
        ('BINY', 'I', None, 'Index for the y coordinate'),
        ('SLOW', 'D', None, 'Parameters for the charging slow component'),
        ('FAST', 'D', None, 'Parameters for the charging fast component')
    ]


def read_charging_map(file_path):
    """ Open a charging map file and read the values for the slow and
    fast component. We get the size of the arrays from the header."""
    logger.info('Reading charging map from %s...', file_path)
    charging_map = fits.open(file_path)['CHRG_MAP']
    nside = charging_map.header['NUM_BINS']
    # Initialize the initial slow and fast values as 2d arrays filled with zeroes
    charging_fast = numpy.full((nside, nside), 0.)
    charging_slow = numpy.full((nside, nside), 0.)
    # We fill the two arrays without doing any assumption on the order of the
    # values in the FITS files, but using explicitly the index from the BINX
    # and BINY columns. This is slower but safer, as it will work even if we
    # change the ordering in the input file.
    data_fast = charging_map.data['FAST']
    data_slow = charging_map.data['SLOW']
    binx = charging_map.data['BINX']
    biny = charging_map.data['BINY']
    for i, (x, y) in enumerate(zip(binx, biny)):
        # Note: indexes of numpy are (row, column), so y goes first
        charging_fast[y, x] = data_fast[i]
        charging_slow[y, x] = data_slow[i]
    logger.info('Done.')
    return charging_fast, charging_slow


def create_charging_map_extension(fast_map, slow_map=None, start_date='N/A',
                                  start_time='N/A', version=1):
    """ Create the CHRG_MAP extension for a charging map file.
    """
    if fast_map.ndim != 2:
        abort('Charging maps must be 2-D arrays.')
    nside = fast_map.shape[0]
    if nside != fast_map.shape[1]:
        abort('Charging map dimensions must be the same on both axes')
    # Note how we use numpy.tile and numpy.repeat to get the right sequence
    # respectively for the rows and the columns
    binx = numpy.tile(numpy.arange(nside), nside)
    biny = numpy.repeat(numpy.arange(nside), nside)
    if slow_map is None:
        slow_map = numpy.full(fast_map.shape, 0.)
    slow = slow_map.flatten()
    fast = fast_map.flatten()
    charging_hdu = xBinTableHDUCharging(data=[binx, biny, slow, fast])
    keywords = {
        'VERSION' : version,
        'CVSD0001' : start_date,
        'CVST0001' : start_time,
        'NUM_BINS' : nside
    }
    for key, value in keywords.items():
        charging_hdu.set_keyword(key, value)
    return charging_hdu


def read_charging_parameters(input_file_path):
    """Open a charging parameters calibration file, read the parameters and
    return them.

    Parameters
    ----------
    input_file_path : string
        the path to the input FITS file storing the charging parameters
    """
    hdu_list = fits.open(input_file_path)
    params = hdu_list['CHRG_PAR'].data[0]
    return params['KC_FAST'], params['TD_FAST'], params['DM_FAST'], \
           params['KC_SLOW'], params['TD_SLOW'], params['DM_SLOW']



class xEnergyFluxCube(xGpdMap3d):

    """Main data structure for the charging-induced gain correction.
    """

    def __init__(self, nside, tedges, zlabel='Mission Elapsed Time [s]',
                 wlabel='Energy flux [keV mm$^{-2}$ s$^{-1}$]'):
        """Overloaded constructor.
        """
        xGpdMap3d.__init__(self, nside, tedges, zlabel=zlabel, wlabel=wlabel)
        # Cache the pixel area, used to normalize the flux at fill() time.
        self.pixel_area = (GPD_PHYSICAL_HALF_SIDE_X / nside) * (GPD_PHYSICAL_HALF_SIDE_Y / nside)
        self.__gain_data = None
        self.charging_map = None

    def xbinning(self):
        """Return the histogram binning on the detx axis.
        """
        return self.binning[0]

    def ybinning(self):
        """Return the histogram binning on the dety axis.
        """
        return self.binning[1]

    def tbinning(self):
        """Return the histogram binning on the time axis (aka z-axis).
        """
        return self.binning[2]

    #pylint: disable=arguments-differ

    def fill(self, detx, dety, time_, energy, deadtime_correction=None, **kwargs):
        """Overloaded fill() method.

        Here we essentially take care of the energy weighting and the
        normalization by the elapsed time and the pixel area. This is not
        completely trivial as the mean energy must be calculated in the
        same time slices used to bin the original histogram, to cope with the
        possibility that the source spectrum---or the event flux on the
        detector, for whatever reason---changes with time.

        Warning
        -------
        The deadtime correction is passed as a single number, i.e., averaged
        over the full data set.
        """
        # Do an initial filling of the three-dimensional histogram with the
        # number of events per space-time cube.
        xGpdMap3d.fill(self, detx, dety, time_, **kwargs)

        # Create a temporary 2-dimensional histogram in time and energy to
        # calculate the average event energy in appropriate time slices.
        tbinning = self.tbinning()
        # Add some padding to the energy binning to deal with the case
        # where we are passing a monochromatic line.
        pad = 0.1
        ebinning = numpy.linspace(energy.min() - pad, energy.max() + pad, 501)
        energy_hist = xHistogram2d(tbinning, ebinning).fill(time_, energy)
        # Calculate the mean energy in each time slice.
        ebin_centers = energy_hist.bin_centers(1)
        weights = energy_hist.content
        energy_sum = weights.sum(axis=1)
        mean_energy = (ebin_centers * weights).sum(axis=1) / energy_sum
        mask = energy_sum > 0.
        mean_energy[~mask] = 0.

        # And now we can calculate the actual energy flux per unit area and
        # overwrite the histogram content. Let's go step by step.
        # First comes the rates in counts per second.
        flux = self.content / numpy.diff(tbinning)

        if deadtime_correction is not None:
            # Next the deadtime correction---note that, in general, the deadtime
            # correction will depend on time, as the counting rate does.
            # Calculate the ovearall *measured* rate on the GPD for each time bin.
            overall_flux = flux.sum(axis=0).sum(axis=0)
            # Translate the average deadtime correction coefficient in the input
            # FITS file into the corresponding average dead time per event.
            dead_time = (1. - deadtime_correction) / overall_flux.mean()
            # Calculate the corrected rate in every time bin using the usual
            # formula (note that overall_rate is an array, here).
            flux /= (1. - overall_flux * dead_time)

        # And, finally, the product with the average energy and the
        # normalization to the pixel area.
        flux *= mean_energy / self.pixel_area

        # And, it goes without saying, update the histogram.
        self.set_content(flux)

    def time_slice(self, tmin=None, tmax=None):
        """Return a time slice of the irradiation history, between generic
        bounds.

        By default this is returning the average energy flux per unit area
        over the entire input sample.

        Note that the times are adjusted to the closest bih edges.
        """
        # Convert the physical times into the corresponding bin edges on the
        # time axis.
        if tmin is None:
            _min = 0
        else:
            _min = numpy.searchsorted(self.tbinning(), tmin)
        if tmax is None:
            _max = -1
        else:
            _max = numpy.searchsorted(self.tbinning(), tmax)
        # Create an empty hyper-histogram to hold the map of the energy flux
        # per unit area.
        hist = xHistogram2d(*self.binning[0:2], *self.labels[0:2], self.labels[3])

        # Do the appropriate averaging of the glorious three-dimensional weights
        # in order to have the appropariate bidimensional slice.
        flux = self.content[:, :, _min:_max].mean(axis=2)
        # Set the contents and return the histogram.
        hist.set_content(flux)
        return hist

    @staticmethod
    def step(flux, dt, dg, k_c, delta_max, tau_d):
        """ Compute a step of the charging model (for either the slow or fast
        process). For a description of the charging parameters see the
        calculate_gain_data() method

        Parameters
        ----------
        flux: float or array
            the value of the input energy flux for this step (ADC counts/mm^2/s)

        dt : float
            the time step (in seconds)

        dg : float or array
            the current gain drop due to this process (fractional)
        """
        if (delta_max == 0) or (k_c == 0):
            return numpy.full(flux.shape, 0.)
        # Add the effect of the charging.
        delta = flux * delta_max / k_c * (1. - dg / delta_max)
        if tau_d != 0:
            # ...and next the discharge.
            delta -= (dg / tau_d)
        return delta * dt

    def calculate_gain_data(self, k_c_fast=DEFAULT_K_C, tau_d_fast=DEFAULT_TAU_D,
        delta_max_fast=DEFAULT_DELTA_MAX, initial_dg_fast=None, k_c_slow=0., tau_d_slow=0.,
        delta_max_slow=0., initial_dg_slow=None):
        """Calculate the expected space-time profile of the gain for a given
        illumination history. We are adopting a charging model with two
        components: one fast and one slow. Both are regulated by the same
        equation, which includes a rate-dependent charging and a constant
        discharging, with the caracteristic time for the two processes (fast
        and slow) separated by roughly an order of magnitude. The system is
        integrated using finite differences. The returned gain profile is
        computed on the times defined by self.tbinning() and has the same
        spatial dimensions as self.content

        Arguments
        ----------
        k_c_fast : float
            the charging constant for the fast process (ADC counts/mm^2)

        tau_d_fast : float
            the discharge time constant for the fast process (seconds)

        delta_max_fast : float
            maximum fractional gain drop for the fast process

        initial_dg_fast : float or array
            initial fraction of gain drop due to the fast process. If an array,
            the value is given for each spatial bin - thus the dimension must
            match the spatial dimensions of self.content. If None, we
            assume initial full discharge for the fast component.

        k_c_slow : float
            as k_c_fast, but for the slow process

        tau_d_slow : float
            as tau_d_fast, but for the slow process

        delta_max_slow : float
            as delta_max_fast, but for the slow process

        initial_dg_slow : float
            as initial_dg_fast, but for the slow process
        """
        # The shape of the gain profile: spatial dimensions as self.content
        # and time axis (the last one) as self.time_grid
        tbinning = self.tbinning()
        gain_shape = self.content.shape[:-1] + (tbinning.shape[0], )
        # initial_dg_fast and initial_dg_slow have the same spatial dimensions
        # of the gain
        if initial_dg_fast is None:
            initial_dg_fast = numpy.full(gain_shape[:-1], 0.)
        if initial_dg_slow is None:
            initial_dg_slow = numpy.full(gain_shape[:-1], 0.)
        # Calculate the initial gain drop across the active surface.
        delta_g_fast = initial_dg_fast * delta_max_fast
        delta_g_slow = initial_dg_slow * delta_max_slow
        # Initialize the gain
        self.__gain_data = numpy.full(gain_shape, 1.)
        # Note: selecting the last axis no matter the rank with numpy '...'
        # syntax
        self.__gain_data[..., 0] = self.__gain_data[..., 0] - delta_g_fast - \
                                   delta_g_slow
        # Cache the time differences for the loop
        delta_tbinning = numpy.diff(tbinning)
        # Time loop. Note the clever way to loop over the last dimentions of the
        # 3d array holding the energy flux, i.e., taking the transpose
        # (by default numpy uses the first axis to iterate).
        # https://stackoverflow.com/questions/1589706/
        # This has no effect for 1-d arrays
        for i, (dt, f) in enumerate(zip(delta_tbinning, self.content.T)):
            # Because of the transpose 'trick' we have to take the transpose
            # of the flux 'f' here, for consistency (otherwise the other
            # dimensions are messed up). See the discussion at
            # https://bitbucket.org/ixpesw/ixpeobssim/issues/389
            # for the need of a transpose of the flux here.
            dg_fast = self.step(f.T, dt, delta_g_fast, k_c_fast, delta_max_fast,
                                tau_d_fast)
            dg_slow = self.step(f.T, dt, delta_g_slow, k_c_slow, delta_max_slow,
                                tau_d_slow)
            delta_g_fast += dg_fast
            delta_g_slow += dg_slow
            self.__gain_data[..., i+1] = self.__gain_data[..., i] - dg_fast - dg_slow
        if delta_max_slow > 0.:
            delta_g_slow /= delta_max_slow
        if delta_max_fast > 0.:
            delta_g_fast /= delta_max_fast
        self.charging_map = [delta_g_fast, delta_g_slow]
        return self.__gain_data

    def gain_map(self, t):
        """Return a two-dimensional histogram with the gain map at a given time.
        """
        hist = xHistogram2d(*self.binning[0:2], *self.labels[0:2], 'Relative gain')
        k = xHistogramBase.bisect(self.tbinning(), t)
        gain = self.__gain_data[:, :, k]
        hist.set_content(gain)
        return hist

    def gain(self, detx, dety, time_, **charging_params):
        """Return the expected corrected gain for a list of events at
        specific detector positions and times.
        """
        if self.__gain_data is None:
            self.calculate_gain_data(**charging_params)
        i = xHistogramBase.bisect(self.xbinning(), detx)
        j = xHistogramBase.bisect(self.ybinning(), dety)
        # Here we specify that the side is 'rigth' to prevent the very first
        # element from being assigned to the last time bin
        k = xHistogramBase.bisect(self.tbinning(), time_, side='right')
        return self.__gain_data[i, j, k + 1]

    def _plot(self):
        """Overloaded plotting method.
        """
        self.time_slice().plot()

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

"""Spurious modulation facilities.
"""

import numbers

import numpy
from scipy.signal import convolve2d
from scipy.stats import norm

from ixpeobssim.core.hist import xHistogram1d, xHistogram2d
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.core.stokes import xModelStokesParameters
from ixpeobssim.core.modeling import xGaussian
from ixpeobssim.core.fitting import fit_histogram
from ixpeobssim.evt.spurmrot import correct_phi_stokes
from ixpeobssim.instrument.gpd import GPD_PHYSICAL_HALF_SIDE_X, NUM_LASER_SWEEPS,\
    LASER_ETCHING_PITCH
from ixpeobssim.irf.base import xResponseBase
from ixpeobssim.irf.modf import xAzimuthalResponseGenerator
from ixpeobssim.irf.ebounds import ENERGY_STEP
from ixpeobssim.utils.matplotlib_ import plt, setup_gca


# pylint: disable=invalid-name, no-member


class xSyntheticSpuriousModulation:

    """Small container class encapsulating a bunch of functions related to the
    simulation of spurious modulation (we're basically using this as a
    namespace).
    """

    @staticmethod
    def amplitude_peak_positions(offset=0.):
        """Return an array of positions (in either the x or y coordinate)
        correponding to the peaks of the spurious modulation, assuming a fixed
        pitch for the laser sweep at the dielectric-drilling stage.

        By default this assumes that the positions are exactly the same in the
        two orthogonal directions, and that the peaks are symmetric with respect
        to the center of the detector itself, although the latter assumption can
        be controlled by the offset parameter (defaulting to zero).
        """
        delta = (NUM_LASER_SWEEPS / 2. - 0.5) * LASER_ETCHING_PITCH
        return numpy.linspace(-delta, delta, NUM_LASER_SWEEPS) + offset

    @staticmethod
    def amplitude_profile(peak_sigma=0.25, min_value=0., max_value=0.1, grid_size=200, offset=0.):
        """Poor man's attempt at building a one-dimensional profile of the
        spurious modulation.

        This is achieved by the sum of a constant and a number of gaussians
        centered on the spurious modulation peaks.
        """
        x = numpy.linspace(-GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_X, grid_size)
        delta = max_value - min_value
        y = numpy.full(x.shape, min_value)
        amplitude = delta * numpy.sqrt(2. * numpy.pi) * peak_sigma
        for peak_pos in xSpuriousModulation.amplitude_peak_positions(offset):
            y += amplitude * norm(peak_pos, peak_sigma).pdf(x)
        fmt = dict(xlabel='x or y [mm]', ylabel='Spurious modulation amplitude')
        return xInterpolatedUnivariateSpline(x, y, **fmt)



class xSpuriousModulationMap:

    """Class describing a spurious modulation map at a given energy.
    """

    NSIDE = 300
    SHAPE = (NSIDE, NSIDE)
    GPD_SIZE = 15.0
    GPD_HALF_SIZE = 0.5 * GPD_SIZE
    CENTER_RADIUS = 58

    def __init__(self, Q, dQ, U, dU):
        """Constructor.
        """
        for a in (Q, dQ, U, dU):
            assert a.shape == self.SHAPE
        self.Q = Q
        self.dQ = dQ
        self.U = U
        self.dU = dU

    @classmethod
    def dummy(cls, Q, U, sigma=0.1):
        """Alternative constructor to create dummy spurious modulation maps with
        given Stokes parameters.
        """
        if isinstance(Q, numbers.Number):
            Q = numpy.full(cls.SHAPE, Q)
        if isinstance(U, numbers.Number):
            U = numpy.full(cls.SHAPE, U)
        Q = Q + numpy.random.normal(0., sigma, size=cls.SHAPE)
        U = U + numpy.random.normal(0., sigma, size=cls.SHAPE)
        dQ = dU = numpy.full(cls.SHAPE, sigma)
        return cls(Q, dQ, U, dU)

    @staticmethod
    def _rebin_array(a, rebin_factor):
        """Custom function to rebin a given map.

        Note this is preserving the original shape of the array, and replacing
        the pixel values in each rebin_factor x rebin_factor patch with the
        average value. (This is so in order not to change the CALDB structure and
        associated code.)

        .. warning::
            There might be a viable, vectorized approach along the lines of
            https://stackoverflow.com/questions/8090229/
            to achieve the same thing, but I don't think this is a bottleneck.
        """
        nx, ny = a.shape
        # A copy of the array is in order, here, not to change the original map.
        a = a.copy()
        # Note the horrible explicit for loop in Python :-(
        for i in range(0, nx, rebin_factor):
            for j in range(0, ny, rebin_factor):
                sel = slice(i, i + rebin_factor), slice(j, j + rebin_factor)
                a[sel] = a[sel].mean()
        return a

    @staticmethod
    def _rebin_uncertainty_array(a, rebin_factor):
        """Specialized facility to propagate the statistical uncertainties in the
        rebin process.
        """
        return numpy.sqrt(xSpuriousModulationMap._rebin_array(a**2., rebin_factor)) / rebin_factor

    def rebin(self, rebin_factor):
        """Return a copy of the map, rebinned by a given factor.
        """
        if rebin_factor <= 1:
            return self
        Q = self._rebin_array(self.Q, rebin_factor)
        dQ = self._rebin_uncertainty_array(self.dQ, rebin_factor)
        U = self._rebin_array(self.U, rebin_factor)
        dU = self._rebin_uncertainty_array(self.dU, rebin_factor)
        return self.__class__(Q, dQ, U, dU)

    @staticmethod
    def _smooth_array_base(a, kernel):
        """Convolve a 2d array with a generic kernel.
        """
        return convolve2d(a, kernel, 'same')

    @staticmethod
    def _smooth_uncertainty_array_base(a, kernel):
        """2d convolution for uncertainties.
        """
        nx, ny = kernel.shape
        return numpy.sqrt(xSpuriousModulationMap._smooth_array_base(a**2., kernel) / (nx * ny))

    @staticmethod
    def _kernel_rect(width, height):
        """Return a rectangular, uniform kernel.
        """
        return numpy.ones((width, height), dtype=float) / (width * height)

    @staticmethod
    def _kernel_square(side):
        """Return a square, uniform kernel.
        """
        return xSpuriousModulationMap._kernel_rect(side, side)

    def _smooth_base(self, kernel):
        """Basic smoothing function.
        """
        if kernel.shape == (1, 1):
            return self
        Q = self._smooth_array_base(self.Q, kernel)
        dQ = self._smooth_uncertainty_array_base(self.dQ, kernel)
        U = self._smooth_array_base(self.U, kernel)
        dU = self._smooth_uncertainty_array_base(self.dU, kernel)
        return self.__class__(Q, dQ, U, dU)

    def smooth_rect(self, width, height):
        """Smooth the map with a rectangular, uniform kernel.
        """
        return self._smooth_base(self._kernel_rect(width, height))

    def smooth_square(self, side):
        """Smooth the map with a square, uniform kernel.
        """
        return self._smooth_base(self._kernel_square(side))

    @classmethod
    def world_to_pixel(cls, detx, dety):
        """Basic conversion from physical to pixel coordinates.

        Note this is implied by the way the spurious modulation maps are created,
        but the thing is not explicitely written in the FITS files.
        """
        if isinstance(detx, numbers.Number):
            detx = numpy.array([detx], float)
        if isinstance(dety, numbers.Number):
            dety = numpy.array([dety], float)
        assert detx.shape == dety.shape
        px = ((detx + cls.GPD_HALF_SIZE) / cls.GPD_SIZE * cls.NSIDE).astype(int)
        py = ((dety + cls.GPD_HALF_SIZE) / cls.GPD_SIZE * cls.NSIDE).astype(int)
        return px, py

    def __call__(self, detx, dety):
        """Overloaded __call__ operator.

        Return the spurious modulation in the form of a 4-element tuple (Q, dQ, U, dU)
        calculated in the physical positions (x, y).
        """
        index = numpy.ravel_multi_index(self.world_to_pixel(detx, dety), self.SHAPE)
        return [a.take(index) for a in (self.Q, self.dQ, self.U, self.dU)]

    def correct_phi(self, phi, detx, dety):
        """Correct an array of azimuthal angles for the spurious modulation.

        .. warning::
            Note this is the old rotation correction that was proven to depend
            on the source parameters, and is therefore fundamentally flawed.
        """
        Q, _, U, dU = self(detx, dety)
        Q += dU**2.
        return correct_phi_stokes(phi, Q, U)

    def ensamble_correction(self, detx, dety):
        """Return the average ensamble correction (Q, U) for a given array of
        detector positions.

        .. warning::
            The error propagation is plain wrong.
        """
        Q, dQ, U, dU = self(detx, dety)
        sqrtn = numpy.sqrt(len(Q))
        return Q.mean(), numpy.sqrt((dQ**2.).mean()) / sqrtn,\
            U.mean(), numpy.sqrt((dU**2.).mean()) / sqrtn

    def rvs_phi(self, detx, dety, Q0=0., U0=0.):
        """Return a random sample of azimuthal angles at given detector positions.
        """
        Q, _, U, _ = self(detx, dety)
        Q += Q0
        U += U0
        modulation = xModelStokesParameters.polarization_degree(Q, U)
        phase = xModelStokesParameters.polarization_angle(Q, U)
        return xAzimuthalResponseGenerator().rvs_phi(modulation, phase)

    @classmethod
    def _fiducial_mask(cls, margin=15):
        """Return a NSIDE x NSIDE boolean mask excluding the borders.
        """
        mask = numpy.zeros((cls.NSIDE, cls.NSIDE), dtype=bool)
        mask[margin:cls.NSIDE - margin, margin:cls.NSIDE - margin] = True
        return mask

    @classmethod
    def _center_mask(cls, radius=None):
        """Return a NSIDE x NSIDE mask selecting a central circle of a given radius.
        """
        if radius is None:
            radius = cls.CENTER_RADIUS
        x, y = numpy.ogrid[:cls.NSIDE, :cls.NSIDE]
        dist = numpy.sqrt((x - cls.NSIDE // 2 + 0.5)**2. + (y - cls.NSIDE // 2 + 0.5)**2.)
        return dist <= radius

    @classmethod
    def _border_mask(cls, radius=None, margin=15):
        """Return a NSIDE x NSIDE mask selecting the complement of the center mask,
        after the application of the fiducial cut.
        """
        return numpy.logical_and(cls._fiducial_mask(margin),\
            numpy.logical_not(cls._center_mask(radius)))

    @classmethod
    def _hist2d(cls, content, zlabel=None):
        """Create a two-dimensional map starting from a generic array.
        """
        binning = numpy.arange(cls.NSIDE + 1)
        fmt = dict(xlabel='X [pixel]', ylabel='Y [pixel]', zlabel=zlabel)
        h = xHistogram2d(binning, binning, **fmt)
        h.set_content(content)
        return h

    def _plot_stokes_params_base(self, data, type_, label='', vmax=0.3):
        """Base function for plotting Stokes maps.
        """
        title = '%s map%s' % (type_, label)
        plt.figure(title)
        self._hist2d(data, zlabel=title).plot(vmin=-vmax, vmax=vmax)
        plt.gca().set_aspect('equal')

    def plot_stokes_params(self, label='', vmax=0.3):
        """Plot the Q and U Stokes parameter maps.
        """
        self._plot_stokes_params_base(self.Q, 'Q', label, vmax)
        self._plot_stokes_params_base(self.U, 'U', label, vmax)

    def plot_error_hist(self, label='', fiducial_mask=True):
        """Plot the one-dimensional histogram of the map errors.

        Note that by default the fiducial mask to remove the map borders is applied.
        """
        binning = numpy.linspace(0., 0.15, 250)
        if fiducial_mask:
            mask = self._fiducial_mask()
        else:
            mask = numpy.full(self.SHAPE, True, dtype=bool)
        title = 'Stokes errors%s' % label
        plt.figure(title)
        dq = self.dQ[mask].flatten()
        du = self.dU[mask].flatten()
        xHistogram1d(binning, xlabel='$\\sigma$').fill(dq).plot(label='Q')
        xHistogram1d(binning, xlabel='$\\sigma$').fill(du).plot(label='U')
        plt.legend()

    def _plot_pulls_base(self, data, type_, label=''):
        """Base function for plotting pulls.
        """
        average = numpy.trunc(numpy.nanmean(data))
        binning = numpy.linspace(average - 7.5, average + 7.5, 100)
        title = '%s pulls%s' % (type_, label)
        plt.figure(title)
        center = self._center_mask()
        border = self._border_mask()
        h = xHistogram1d(binning, xlabel='Pulls [$\\sigma$]').fill(data[border].flatten())
        h.plot(label='Border')
        border_model = fit_histogram(xGaussian(), h)
        border_model.plot()
        border_model.stat_box(position='upper left')
        h = xHistogram1d(binning, xlabel='Pulls [$\\sigma$]').fill(data[center].flatten())
        h.plot(label='Center')
        center_model = fit_histogram(xGaussian(), h)
        center_model.plot()
        center_model.stat_box(position='lower left')
        plt.legend()
        return border_model, center_model

    def plot_pulls(self, label=''):
        """Plot the pulls.
        """
        qbm, qcm = self._plot_pulls_base(self.Q / self.dQ, 'Q', label)
        ubm, ucm = self._plot_pulls_base(self.U / self.dU, 'U', label)
        return qbm, qcm, ubm, ucm

    def plot_correlation(self, label=''):
        """Plot the correlation between Q and U.
        """
        title = 'Correlation%s' % label
        plt.figure(title)
        binning = numpy.linspace(-0.5, 0.5, 50)
        mask = self._fiducial_mask()
        q = self.Q[mask].flatten()
        u = self.U[mask].flatten()
        h = xHistogram2d(binning, binning, xlabel='Q', ylabel='U').fill(q, u)
        h.plot()
        plt.gca().set_aspect('equal')



class xSpuriousModulation(xResponseBase):

    """Class describing a set of spurious modulation maps.
    """

    NUM_LAYERS = 6

    def __init__(self, file_path, extension):
        """Overloaded constructor
        """
        xResponseBase.__init__(self, file_path, 'fits')
        data = self.hdu_list['MODULATION'].data
        Q = self._reshape(data['DETQ_SM'])
        dQ = self._reshape(data['D_DETQ_SM'])
        U = self._reshape(data['DETU_SM'])
        dU = self._reshape(data['D_DETU_SM'])
        self._map_list = []
        for i in range(self.NUM_LAYERS):
            map_ = xSpuriousModulationMap(Q[:, :, i], dQ[:, :, i], U[:, :, i], dU[:, :, i])
            self._map_list.append(map_)
        self._pi_grid = data['PI'][0]
        self._energy_grid = ENERGY_STEP * self._pi_grid

    @classmethod
    def _reshape(cls, a):
        """Reshape a generic array, corresponding to a column in the input FITS
        file, to an actual binned map in (DETX, DETY, PI) space.
        """
        return a.reshape((*xSpuriousModulationMap.SHAPE, cls.NUM_LAYERS))

    def map_(self, layer):
        """Return the spurious modulation map at a given energy layer.
        """
        return self._map_list[layer]

    def _plot_label(self, layer):
        """Return the proper label to be appended to the plot title for a given layer.
        """
        return ' @ %.2f keV' % self._energy_grid[layer]

    def plot_stokes_map(self, layer, vmax=0.3):
        """Plot the underlying Q and U maps in a given energy layer.
        """
        self.map_(layer).plot_stokes_params(self._plot_label(layer), vmax)

    def plot_all_stokes_maps(self, vmax=0.3):
        """Plot all the Q and U maps.
        """
        for layer in range(self.NUM_LAYERS):
            self.plot_stokes_map(layer, vmax)

    def plot_errors(self, layer):
        """Plot the Q and U statistical errors for a given layer.
        """
        self.map_(layer).plot_error_hist(self._plot_label(layer))

    def plot_all_errors(self):
        """Plot the Q and U statistical errors for all the layers.
        """
        for layer in range(self.NUM_LAYERS):
            self.plot_errors(layer)

    def plot_pulls(self, layer):
        """Plot the pulls for a given energy layer.
        """
        return self.map_(layer).plot_pulls(self._plot_label(layer))

    def plot_all_pulls(self):
        """Plot all the Q and U pulls.
        """
        pull_sigma_border_q = []
        pull_sigma_center_q = []
        pull_sigma_border_u = []
        pull_sigma_center_u = []
        for layer in range(self.NUM_LAYERS):
            border_model_q, center_model_q, border_model_u, center_model_u =\
                self.plot_pulls(layer)
            pull_sigma_border_q.append(border_model_q.Sigma)
            pull_sigma_center_q.append(center_model_q.Sigma)
            pull_sigma_border_u.append(border_model_u.Sigma)
            pull_sigma_center_u.append(center_model_u.Sigma)
        plt.figure('Pulls summary')
        plt.plot(self._energy_grid, pull_sigma_center_q, 'o', label='Q center')
        plt.plot(self._energy_grid, pull_sigma_center_u, 'o', label='U center')
        plt.plot(self._energy_grid, pull_sigma_border_q, 'o', label='Q border')
        plt.plot(self._energy_grid, pull_sigma_border_u, 'o', label='U border')
        setup_gca(xlabel='Energy [keV]', ylabel='$\\sigma$ pulls', grids=True, legend=True, ymin=0.)
        for y, label in ((1., 'Random noise'), (2.**0.5, '100% "relative error"'),
            (11.**0.5, '10% "relative error"')):
            plt.axhline(y, ls='dashed', color='black')
            plt.text(5., y + 0.025, label, ha='center')

    def plot_correlation(self, layer):
        """Plot the correlation between Q and U in a given energy layer.
        """
        return self.map_(layer).plot_correlation(self._plot_label(layer))

    def plot_all_correlations(self):
        """Plot the correlation between Q and U in all the energy layers.
        """
        for layer in range(self.NUM_LAYERS):
            self.plot_correlation(layer)



if __name__ == '__main__':
    m = xSpuriousModulationMap.dummy(0., 0., 0.1)
    m.plot_stokes_params(' original')
    m.plot_pulls(' original')
    m1 = m.smooth_square(5)
    m1.plot_stokes_params(' smooth')
    m1.plot_pulls(' smooth')
    m2 = m.smooth_rect(3, 10)
    m2.plot_stokes_params(' smooth 2')
    m2.plot_pulls(' smooth 2')
    plt.show()

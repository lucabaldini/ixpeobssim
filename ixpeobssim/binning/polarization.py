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

"""Binned data structures pertaining to the (spectro-)polarimetric analysis.
"""

from __future__ import print_function, division

import os
import numbers

import astropy
from astropy.io import fits
import numpy
import scipy

from ixpeobssim.binning.base import xEventBinningBase, xBinnedFileBase
from ixpeobssim.binning.fmt import xBinTableHDUPHA1, xBinTableHDUPCUBE, xBinTableHDUEBOUNDS
from ixpeobssim.core.fitsio import xFITSImageBase
from ixpeobssim.core.hist import xScatterPlot
from ixpeobssim.core.stokes import xModelStokesParameters
from ixpeobssim.evt.align import align_stokes_parameters
from ixpeobssim.evt.kislat2015 import xStokesAnalysis
from ixpeobssim.irf import load_modf
from ixpeobssim.irf.caldb import irf_file_path
from ixpeobssim.irf.ebounds import NUM_CHANNELS
from ixpeobssim.utils.astro import region_compound
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, nlog_errorbars,\
    plot_arrows, plot_ellipse, setup_gca_stokes
from ixpeobssim.utils.units_ import arcmin_to_arcsec


# pylint: disable=invalid-name, no-member, attribute-defined-outside-init
# pylint: disable=too-many-locals, too-many-instance-attributes


class xEventBinningPHA1Base(xEventBinningBase):

    """Base class for PHA1 binning.
    """

    ANCRFILE = None
    XFLT0001 = None
    SUPPORTED_KWARGS = ['mc', 'weights', 'weightcol', 'grayfilter']

    def __init__(self, file_path, **kwargs):
        """Overloaded constructor.
        """
        xEventBinningBase.__init__(self, file_path, **kwargs)
        if self.irf_name is None:
            abort('Please set --irfname command-line switch explicitely')
        self.binning = numpy.linspace(0, NUM_CHANNELS, NUM_CHANNELS + 1) - 0.5
        self.channels = numpy.arange(NUM_CHANNELS)
        self.livetime = self.event_file.livetime()

    def _binned_data(self, stokes_weights=None):
        """Base binning algorithm to be specialized by derived classes.

        Args
        ----
        stokes_weights : array_like, optional
            The weights to be assigned to the histogram entries to distinguish
            between the PHA1 (None), PHA1U (2. * sin(2. * phi)) and
            PHA1Q (2. * cos(2. * phi)) binning mode. Note these should not
            be confused with the actual event weights in an ensemble weighted
            analysis---the latter are automaically picked up in the function body
            depending on the command-line arguments passed to the app.

        Note that, in terms of assigning errors to the output spectra when the
        weights are not None, we use the standard formula for weighted
        histograms, where the errors are just the square root of the sum of the
        weights squared, see e.g.
        www.hep.uiuc.edu/e687/memos/weight_err.ps
        """
        # Retrieve the pi.
        pi = self.event_file.pi_data(self.get('mc'))
        # Stokes weights, i.e., I, Q or U?
        if stokes_weights is None:
            stokes_weights = numpy.full(pi.shape, 1.)
        # Event weights.
        event_weights = self.weight_data()
        # Full weights.
        weights = stokes_weights * event_weights
        # Accumulate the necessary histograms.
        w, _ = numpy.histogram(pi, bins=self.binning, weights=weights)
        w2, _ = numpy.histogram(pi, bins=self.binning, weights=weights**2)
        ew, _ = numpy.histogram(pi, bins=self.binning, weights=event_weights)
        ew2, _ = numpy.histogram(pi, bins=self.binning, weights=event_weights**2.)
        # Calculate the "normalization" for our weight prescription.
        # Note the protection against zero-division errors, see
        # https://bitbucket.org/ixpesw/ixpeobssim/issues/470/runtime-warning-in-xpbin
        N = numpy.zeros(ew.shape)
        mask = ew2 != 0
        N[mask] = ew[mask] / ew2[mask] / self.livetime
        # And we're good to go :-)
        rate = N * w
        error = N * numpy.sqrt(w2)
        return self.channels, rate, error

    def _normalized_binned_data(self, weights):
        """Normalized version of the binning routine for Q/I and U/I normalized
        Stokes parameters.

        Note that some care must be taken in performing the division and
        propagating the uncertainties in order to avoid zero-division runtime
        errors.

        For the uncertainties we carefully avoid propagating the relative errors,
        to avoid problems with points where uq == 0, duq != 0 and duq/uq is nan.

        For the bins where i == 0 we currently explicitly set both the value and
        the uncertainty of the normalized Stokes parameter to 0, whith the
        intent of doing the right thing when we rebin spectra.

        Warning
        -------
        We might want to put some more thought into what the error should be
        when i == 0, because in this case, if we don't rebin, we have potential
        problems downstream with the chisquare calculation (although one might
        argue that in this case the chisquare is not the right statistics to
        start with.)
        """
        # Retrieve the releval un-normalized spectra.
        _, uq, duq = self._binned_data(weights)
        _, i, di = self._binned_data(None)
        # Initialize the normalized rate and associated errors.
        rate = numpy.full(i.shape, 0.)
        error = numpy.full(i.shape, 0.)
        # Mask to avoid runtime errors due to zero division downstream.
        mask = i > 0.
        # Restric all the calculation to the bins with at least 1 I count.
        uq = uq[mask]
        duq = duq[mask]
        i = i[mask]
        di = di[mask]
        # Overerwrite with the proper values the portion of the output array
        # where i > 0.
        #
        # Note that we propagate the absolute error as, if we were to
        # perform the error propagation through the relative errors, we
        # could be bitten by occasional data points where uq == 0, duq/uq is nan,
        # and yet the point itself is perfectly legitimate.
        rate[mask] = uq / i
        error[mask] = numpy.sqrt((uq * di)**2. + (i * duq)**2.) / (i**2.)
        return self.channels, rate, error

    def binned_data(self):
        """Method to be overloaded by derived classes.
        """
        raise NotImplementedError

    def bin_(self):
        """Overloaded method.
        """
        primary_hdu = self.build_primary_hdu()
        data = self.binned_data()
        spec_hdu = xBinTableHDUPHA1(data)
        spec_hdu.setup_header(self.event_file.primary_keywords())
        du_id = self.event_file.du_id()
        keywords = [
            ('EXPOSURE', self.event_file.livetime()),
            ('RESPFILE', irf_file_path(self.irf_name, du_id, 'rmf')),
            ('ANCRFILE', irf_file_path(self.irf_name, du_id, self.ANCRFILE,
                gray_filter=self.get('grayfilter'))),
            ('XFLT0001', self.XFLT0001)
            ]
        spec_hdu.setup_header(keywords)
        hdu_list = fits.HDUList([primary_hdu, spec_hdu])
        self.write_output_file(hdu_list)



class xEventBinningPHA1(xEventBinningPHA1Base):

    """Original algorithm for PHA1 files.
    """

    INTENT = 'count spectrum'
    ANCRFILE = 'arf'
    XFLT0001 = 'Stokes:0'

    def binned_data(self):
        """Overloaded binning algorithm.
        """
        return xEventBinningPHA1Base._binned_data(self, None)



class xEventBinningPHA1Q(xEventBinningPHA1Base):

    """Class for PHA1 binning.
    """

    INTENT = 'Stokes Q spectrum'
    ANCRFILE = 'mrf'
    XFLT0001 = 'Stokes:1'

    def binned_data(self):
        """Overloaded binning algorithm.
        """
        return xEventBinningPHA1Base._binned_data(self, self.event_file.q_data())



class xEventBinningPHA1U(xEventBinningPHA1Base):

    """Subclass for creating Stokes U spectra.
    """

    INTENT = 'Stokes U spectrum'
    ANCRFILE = 'mrf'
    XFLT0001 = 'Stokes:2'

    def binned_data(self):
        """Overloaded binning algorithm.
        """
        return xEventBinningPHA1Base._binned_data(self, self.event_file.u_data())



class xEventBinningPHA1QN(xEventBinningPHA1Q):

    """Subclass for creating normalized Stokes Q spectra.
    """

    INTENT = 'normalized Stokes Q spectrum (Q/I)'
    ANCRFILE = 'modf'

    def binned_data(self):
        """Overloaded binning algorithm.
        """
        return self._normalized_binned_data(self.event_file.q_data())



class xEventBinningPHA1UN(xEventBinningPHA1U):

    """Subclass for creating normalized Stokes U spectra.
    """

    INTENT = 'normalized Stokes U spectrum (U/I)'
    ANCRFILE = 'modf'

    def binned_data(self):
        """Overloaded binning algorithm.
        """
        return self._normalized_binned_data(self.event_file.u_data())



class xBinnedCountSpectrum(xBinnedFileBase):

    """Binned count spectrum.
    """

    def _read_data(self):
        """Overloaded method.
        """
        # pylint: disable=no-member
        self.spectrum_header = self.hdu_list['SPECTRUM'].header
        self._read_binary_table_data(xBinTableHDUPHA1.NAME)

    def __iadd__(self, other):
        """ Overloaded method for PHA1 binned file addition.
        """
        # pylint: disable=no-member, attribute-defined-outside-init
        self._check_iadd(other, ('RATE', 'STAT_ERR'), ('CHANNEL', ))
        self.RATE += other.RATE
        self.STAT_ERR = numpy.sqrt(self.STAT_ERR**2. + other.STAT_ERR**2.)
        return self

    def respfile(self):
        """Return the value of the RESPFILE keyword in the SPECTRUM extension
        header, stripped from the first part of the pathe (i.e., the file name
        only).
        """
        return os.path.basename(self.spectrum_header['RESPFILE'])

    def plot(self, **kwargs):
        """Overloaded plot method.
        """
        # pylint: disable=no-member, arguments-differ
        ylabel = kwargs.pop('ylabel', 'Rate [Hz]')
        grids = kwargs.pop('grids', True)
        logy = kwargs.pop('logy', True)
        kwargs.setdefault('fmt', 'o')
        nlog_errorbars(self.CHANNEL, self.RATE, self.STAT_ERR, **kwargs)
        setup_gca(xlabel='Channel', xmax=self.CHANNEL.max(), ylabel=ylabel, grids=grids, logy=logy)


class xEventBinningPCUBE(xEventBinningBase):

    """Class for PCUBE binning.

    .. versionadded:: 12.0.0
    """

    INTENT = 'polarization cube (in energy layers)'
    SUPPORTED_KWARGS = ['mc', 'acceptcorr', 'weights', 'weightcol', 'grayfilter'] + \
        xEventBinningBase._energy_binning_kwargs()

    def bin_(self):
        """Overloaded method.
        """
        if self.irf_name is None:
            abort('Please set --irfname command-line switch explicitely')
        energy = self.event_file.energy_data(self.get('mc'))
        ebinning = self.make_energy_binning(energy, **self.kwargs)
        q, u = self.event_file.stokes_data()
        modf = load_modf(self.irf_name, self.event_file.du_id())
        aeff = self.load_aeff_for_polarization_analysis()
        analysis = xStokesAnalysis(q, u, energy, modf, aeff, self.event_file.livetime(),
            self.weight_data(), self.get('acceptcorr'))
        table = analysis.polarization_table(ebinning)
        data = [table[key] for key in xBinTableHDUPCUBE.spec_names()]
        primary_hdu = self.build_primary_hdu()
        bin_table_hdu = xBinTableHDUPCUBE(data)
        bin_table_hdu.setup_header(self.event_file.primary_keywords())
        gti_hdu = self.event_file.hdu_list['GTI']
        hdu_list = fits.HDUList([primary_hdu, bin_table_hdu, gti_hdu])
        self.write_output_file(hdu_list)



class xBinnedPolarizationCube(xBinnedFileBase):

    """Read-mode interface to a PCUBE FITS file.

    .. versionadded:: 12.0.0
    """

    def _read_data(self):
        """Overloaded method.
        """
        self._read_binary_table_data(xBinTableHDUPCUBE.NAME)

    def backscal(self):
        """Return the value of the BACKSCAL header keyword, if present.
        """
        try:
            return self.primary_header['BACKSCAL']
        except KeyError:
            logger.warning('Polarization cube has no BACKSCAL header keyword set')
            return None

    def __check_compat(self, other):
        """Check the basic polarization cube data structure before attempting
        to do operations with other polarization cubes.
        """
        same_shape = xBinTableHDUPCUBE.POL_COL_NAMES
        same_values = ('ENERG_LO', 'ENERG_HI')
        self._check_iadd(other, same_shape, same_values)

    def __recalculate_derived(self):
        """Recalculate all the derived quantities after arithmetic operations.
        """
        # Recalculate the normalized Stokes parameters, and propagate the errors.
        args = self.I, self.Q, self.U, self.MU, self.W2
        self.QN, self.UN, self.I_ERR, self.Q_ERR, self.U_ERR, self.QN_ERR, self.UN_ERR, \
            self.QUN_COV, self.P_VALUE, self.CONFID, self.SIGNIF = \
                xStokesAnalysis.calculate_stokes_errors(*args)
        # Update the MDP.
        self.MDP_99 = xStokesAnalysis.calculate_mdp99(self.MU, self.I, self.W2)
        self.N_EFF, self.FRAC_W = xStokesAnalysis.calculate_n_eff(self.COUNTS, self.I, self.W2)
        # Recalculate the polarization degree and angle.
        self.PD, self.PD_ERR, self.PA, self.PA_ERR = \
            xStokesAnalysis.calculate_polarization(*args, degrees=True)

    def __iadd__(self, other):
        """Overloaded method for PCUBE binned data addition.
        """
        self.__check_compat(other)
        self.E_MEAN = self._weighted_average(other, 'E_MEAN', 'I')
        self.MU = self._weighted_average(other, 'MU', 'I')
        self.COUNTS += other.COUNTS
        self.W2 += other.W2
        self.I += other.I
        self.Q += other.Q
        self.U += other.U
        self.__recalculate_derived()
        return self

    def __isub__(self, other):
        """Overloaded method for polarization cube subtraction.
        """
        self.__check_compat(other)
        self.E_MEAN = self._weighted_average(other, 'E_MEAN', 'I', invert_w2=True)
        self.MU = self._weighted_average(other, 'MU', 'I', invert_w2=True)
        self.COUNTS -= other.COUNTS
        # Note W2 must always be added.
        self.W2 += other.W2
        self.I -= other.I
        self.Q -= other.Q
        self.U -= other.U
        # Propagate the uncertainties in the polarization cube difference.
        # This was added in response to https://bitbucket.org/ixpesw/ixpeobssim/issues/614
        # and should be considered as a short-term fix, waiting for a full
        # refactoring of the polarization cube arithmetic.
        self.I_ERR = numpy.sqrt(self.I_ERR**2. + other.I_ERR**2.)
        self.Q_ERR = numpy.sqrt(self.Q_ERR**2. + other.Q_ERR**2.)
        self.U_ERR = numpy.sqrt(self.U_ERR**2. + other.U_ERR**2.)
        args = self.I, self.Q, self.U, self.I_ERR, self.Q_ERR, self.U_ERR
        self.QN, self.UN, self.QN_ERR, self.UN_ERR = \
            xStokesAnalysis.normalized_stokes_parameters(*args)
        self.QUN_COV = xStokesAnalysis.stokes_covariance(self.I, self.QN, self.UN, self.W2)
        self.P_VALUE, self.CONFID, self.SIGNIF = \
            xStokesAnalysis.significance(self.Q, self.U, self.Q_ERR, self.U_ERR)
        self.MDP_99 = xStokesAnalysis.calculate_mdp99(self.MU, self.I, self.W2)
        self.N_EFF, self.FRAC_W = xStokesAnalysis.calculate_n_eff(self.COUNTS, self.I, self.W2)
        self.PD, self.PD_ERR, self.PA, self.PA_ERR = \
            xStokesAnalysis.calculate_polarization_sub(*args, degrees=True)
        return self

    def __imul__(self, other):
        """Overloaded method for polarization cube multiplication by a scalar.
        """
        assert isinstance(other, numbers.Number)
        # Need to be careful, here, as the counts are supposed to be integer.
        counts = (self.COUNTS * other + 0.5).astype(int)
        self.COUNTS = counts
        # Note the square in W2!
        self.W2 *= other**2
        self.I *= other
        self.Q *= other
        self.U *= other
        self.__recalculate_derived()
        return self

    def polarization(self):
        """Return the polarization information in the form of a four element
        tuple (PD, PD_err, PA, PA_err) of arrays.
        """
        return self.PD, self.PD_ERR, self.PA, self.PA_ERR

    def _energy_scatter_plot(self, y, dy, ylabel=None):
        """Make a scatter plot of a given quanity vs. energy.

        Basic function for plot_polarization_degree() and plot_polarization_angle().
        """
        x = self.E_MEAN
        dx = [self.E_MEAN - self.ENERG_LO, self.ENERG_HI - self.E_MEAN]
        return xScatterPlot(x, y, dy, dx, xlabel='Energy [keV]', ylabel=ylabel)

    def plot_polarization_degree(self, min_num_sigmas=2., **kwargs):
        """Plot the polarization degree as a function of energy.

        .. warning::

           This is deprecated in favor of the standard plots of the Stokes parameters.
        """
        ylabel = 'Polarization degree'
        limits = (self.PD / self.PD_ERR) < min_num_sigmas
        PD = self.PD + limits * min_num_sigmas * self.PD_ERR
        scatter = self._energy_scatter_plot(PD, self.PD_ERR, ylabel)
        scatter.plot(uplims=limits, **kwargs)

    def plot_polarization_angle(self, **kwargs):
        """Plot the polarization angle as a function of energy.

        .. warning::

           This is deprecated in favor of the standard plots of the Stokes parameters.
        """
        ylabel = 'Polarization angle [$^\\circ$]'
        scatter = self._energy_scatter_plot(self.PA, self.PA_ERR, ylabel)
        scatter.plot(**kwargs)

    def plot(self, **kwargs):
        """Default plotting for the polarization cubes.

        .. warning::

           This should be considered in evolution, and the API might change.

        This is plotting the polarization cube in normalized Stokes parameters
        space, with custom grids to help reading out the polarization degree and
        angle.

        A brief explanation of the supporte keywords arduments.

        Arguments
        ---------
        side : float, optional
            The absolute value of the maximum Q/I and U/I to be displayed
            (note the aspect ratio of the plot is set to 'equal', and the
            plot is assumed to be squared). By default this is driven by the
            maximum value of the polarization degree over the polarization cube.

        pd_grid : array_like
            The polarization degree radial grid in correspondence of which the
            gridding circumpherences are plotted to guide the eye.

        pd_grid_label_angle : float
            The angle at which the text labels for the PD grids are plotted.

        pa_grid_step : float
            The step of the tangential grid in the polarization angle.

        sigma_levels : array_like
            The sigma levels at which the ellipses repesenting the normalized
            Stokes parameter contours are plotted.

        sigma_ls : tuple of the same size as sigma_levels
            The line styles corresponding to sigma levels---by default the
            innermost ellipse is solid and the others are dashed.

        colors : array_like or str, optional
            The colors for the elliptical contours (by deafult the proper color
            for the specific DU in the polarization cube is picked).

        label : str, optional
            Optional label to be displayed in the legend---this is only used if
            the marker boolean flag is True and only associated to the first
            layer in the cube.

        marker : bool
            If True, a marker is plotted at the center of the elliptical contours.

        marker_size : float
            The size for the optional markers.

        annotate_energies : bool
            If true, the energy layers are annotated with nice arrows and text
            labels (note the positioning is still fragile).

        setup_axes : bool
            Call the underlying setup_gca_stokes() hook at the end. (This flag
            is provided to allow disingaging multiple grid plotting in case one
            wants to overlay multiple polarization cubes.)
        """
        # Cache all the relevant command-line arguments.
        side = kwargs.get('side', None)
        pd_grid = kwargs.get('pd_grid')
        pd_grid_label_angle = kwargs.get('pd_grid_label_angle', 45.)
        pa_grid_step = kwargs.get('pa_grid_step', 30.)
        sigma_levels = kwargs.get('sigma_levels', (1., 2., 3.))
        sigma_ls = kwargs.get('sigma_ls', None)
        colors = kwargs.get('colors', None)
        label = kwargs.get('label', None)
        marker = kwargs.get('marker', True)
        marker_size = kwargs.get('marker_size', 4)
        annotate_energies = kwargs.get('annotate_energies', True)
        setup_axes = kwargs.get('setup_axes', True)

        # Setup the default side---this is just the maximum of the polarization
        # degree across the cube plus twice the associated sigma.
        if side is None:
            index = numpy.argmax(self.PD)
            side = min(self.PD[index] + 2. * self.PD_ERR[index], 1.)
        # Setup the radial grid for the polarization angle. In order to have a
        # sensible default, here, we judge based on the plot side and try and
        # setup a sensible regular grid.
        if pd_grid is None:
            if side <= 0.05:
                step = 0.01
            elif side <= 0.1:
                step = 0.02
            elif side <= 0.25:
                step = 0.05
            else:
                step = 0.1
            pd_grid = numpy.arange(step, side, step)
        # If the color is not set, default to the color for the DU.
        if colors is None:
            colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        elif isinstance(colors, str):
            colors = [colors] * len(self.QN)
        # Default line style is solid for the first level, and dashed for the others.
        if sigma_ls is None:
            sigma_ls = ['solid'] + ['dashed'] * (len(sigma_levels) - 1)

        # Loop over the actual data.
        args = self.QN, self.UN, self.QN_ERR, self.UN_ERR, self.ENERG_LO, self.ENERG_HI, colors
        for i, (q, u, dq, du, emin, emax, color) in enumerate(zip(*args)):
            # Plot all the contours relative to the desired sigma levels.
            for scale, ls in zip(sigma_levels, sigma_ls):
                plot_ellipse((q, u), 2. * scale * dq, 2. * scale * du, color=color,
                    ls=ls, zorder=10)
            # Plot the markers.
            if marker:
                _kwargs = dict(ms=marker_size, color=color, zorder=10)
                # For the first layer, we also set the label, so that we can
                # effectively use a legend downstream.
                if i == 0:
                    plt.plot(q, u, 'o', label=label, **_kwargs)
                else:
                    plt.plot(q, u, 'o', **_kwargs)
            # Annotations for the energy layers.
            if annotate_energies:
                label = '%.2f-%.2f keV' % (emin, emax)
                if q > 0:
                    xtext = q - 0.05 + 0.02 * i
                else:
                    xtext = q + 0.05 - 0.02 * i
                if u > 0:
                    ytext = u - 0.05 - 0.02 * i
                else:
                    ytext = u + 0.05 + 0.02 * i
                plt.gca().annotate(label, xy=(q, u), xycoords='data', xytext=(xtext, ytext),
                    textcoords='data', size='small', ha='center', zorder=10, color=color,
                    arrowprops=dict(arrowstyle="->", color=color, lw=0.75,
                    shrinkA=5, shrinkB=1, patchA=None, patchB=None,
                    connectionstyle='angle3,angleA=90,angleB=0'))
        # Finally, set up the axes.
        if setup_axes:
            setup_gca_stokes(side, pd_grid, pd_grid_label_angle, pa_grid_step)

    def as_table(self):
        """Return the polarization cube as an astropy table.
        """
        col_names = ['Quantity']
        col_names += ['%.2f--%.2f keV' % (emin, emax) for emin, emax in \
            zip(self.ENERG_LO, self.ENERG_HI)]
        col_types = [str] + [float] * (len(col_names) - 1)
        table = astropy.table.Table(names=col_names, dtype=col_types)
        for col_name, *_ in xBinTableHDUPCUBE.DATA_SPECS[2:]:
            table.add_row([col_name] + list(self.__getattr__(col_name)))
        return table

    def __str__(self):
        """String formatting.
        """
        return '%s content:\n%s' % (self.__class__.__name__, self.as_table())



class xEventBinningMDPMAPCUBE(xEventBinningBase):

    """Class for MDPMAPCUBE binning.
    """

    INTENT = 'MDP map cube in sky coordinates and energy layers'
    SUPPORTED_KWARGS = ['mc', 'acceptcorr', 'weights', 'weightcol', 'grayfilter'] + \
        xEventBinningBase._image_wcs_kwargs() + xEventBinningBase._energy_binning_kwargs()
    EXTENSIONS_NAMES = xBinTableHDUPCUBE.MDP_COL_NAMES

    def process_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.process_kwargs(self)
        self.process_image_ref_kwargs()

    @staticmethod
    def _binned_data(stokes_analysis, *args, **kwargs):
        """Retrieve the necessary data from the underlying xStokesAnalysis object.
        """
        return stokes_analysis.mdp_map_cube(*args, **kwargs)

    def bin_(self):
        """Overloaded method.
        """
        if self.irf_name is None:
            abort('Please set --irfname command-line switch explicitely')
        # Retrieve the relevant columns in the input event file.
        ra, dec = self.event_file.sky_position_data(self.get('mc'))
        energy = self.event_file.energy_data(self.get('mc'))
        ebinning = self.make_energy_binning(energy, **self.kwargs)
        q, u = self.event_file.stokes_data()
        # Load the necessary response functions.
        modf = load_modf(self.irf_name, self.event_file.du_id())
        aeff = self.load_aeff_for_polarization_analysis()
        analysis = xStokesAnalysis(q, u, energy, modf, aeff, self.event_file.livetime(),
            self.weight_data(), self.get('acceptcorr'))
        # Prepare the arrays and binning for mapping the polarization.
        wcs_ = self._build_image_wcs(**self.kwargs)
        x, y, xybinning = self._pixelize_skycoords(ra, dec, wcs_)
        # Note the energy axis is up front in order to be able to view the
        # images with fv at a later stage.
        binning = (ebinning, *xybinning)
        binned_data = self._binned_data(analysis, x, y, binning)
        # Prepare the output HDU list.
        header = wcs_.to_header()
        hdu_list = fits.HDUList()
        hdu = self.build_primary_hdu(header=header)
        hdu_list.append(hdu)
        for ext_name in self.EXTENSIONS_NAMES:
            data = binned_data[ext_name]
            hdu_list.append(fits.ImageHDU(data, header=header, name=ext_name))
        logger.info('Creating EBOUNDS...')
        ebounds = xBinTableHDUEBOUNDS((ebinning[:-1], ebinning[1:]))
        ebounds.setup_header(self.event_file.primary_keywords())
        hdu_list.append(ebounds)
        self.write_output_file(hdu_list)



class xBinnedMDPMapCube(xBinnedFileBase):

    """Read-mode interface to a MDPMAPCUBE FITS file.
    """

    EXTENSIONS_NAMES = xEventBinningMDPMAPCUBE.EXTENSIONS_NAMES

    def _read_data(self):
        """Overloaded method.
        """
        self._img_header = self.hdu_list['I'].header
        self.wcs = astropy.wcs.WCS(self._img_header)
        for ext_name in self.EXTENSIONS_NAMES:
            self._read_image_data(ext_name)
        self._read_binary_table_data(xBinTableHDUEBOUNDS.NAME)

    def __iadd__(self, other):
        """Overloaded method for binned data addition.
        """
        same_shape = xBinTableHDUPCUBE.MDP_COL_NAMES
        same_values = ('ENERG_LO', 'ENERG_HI')
        self._check_iadd(other, same_shape, same_values)
        self.E_MEAN = self._weighted_average(other, 'E_MEAN', 'I')
        self.COUNTS += other.COUNTS
        self.MU = self._weighted_average(other, 'MU', 'I')
        self.W2 += other.W2
        self.I += other.I
        self.MDP_99 = xStokesAnalysis.calculate_mdp99(self.MU, self.I, self.W2)
        self.N_EFF, self.FRAC_W = xStokesAnalysis.calculate_n_eff(self.COUNTS, self.I, self.W2)
        return self

    def map_shape(self):
        """Return the shape of the underlying sky-maps.

        Mind the underlying arrays are all 3-dimensional cubes with the energy
        binning as the first dimension---so it's the last two that we care about.
        (And, for completeness: since the arrays are all the same, we use I
        for convenience.)
        """
        return self.I.shape[1:]

    def pixel_size(self):
        """Return the pixel size of the underlying spatial map.
        """
        cdelt1, cdelt2, _ = self.wcs.wcs.get_cdelt()
        assert abs(cdelt1) == abs(cdelt2)
        return abs(cdelt1)


    def ds9_region_mask(self, region_list, region_slice=None):
        """Return the (spatial) array map corresponding to a given ds9 region list.

        If there is more than one region in the region list, by default the
        output mask corresponds to the logical or of all the regions.
        The region_slice optional argument allows to restrict the mask to a
        subset of the regions.

        Arguments
        ---------
        region_list :
            The region list defining the mask

        region_slice : int or slice (optional)
            An optional slice designator to select a subset of the region list.
        """
        if region_slice is not None:
            if isinstance(region_slice, int):
                region_slice = slice(region_slice, region_slice + 1)
            region_list = region_list[region_slice]
        region = region_compound(region_list)[0].to_pixel(self.wcs)
        return region.to_mask().to_image(self.map_shape()).astype(bool)

    @staticmethod
    def _sum_image_pixels(data, spatial_mask, energy_layer=0):
        """Sum one of the underlying image extentions over a given spatial
        mask and energy slice.
        """
        return data[energy_layer][spatial_mask].sum()

    def sum_pixels(self, spatial_mask, energy_layer=0):
        """Sum the relevant quantities over a given spatial mask and energy slice.
        """
        args = spatial_mask, energy_layer
        I = self._sum_image_pixels(self.I, *args)
        counts = self._sum_image_pixels(self.COUNTS, *args)
        mu = self._sum_image_pixels(self.MU * self.I, *args)
        if I > 0.:
            mu /= I
        W2 = self._sum_image_pixels(self.W2, *args)
        mdp = xStokesAnalysis.calculate_mdp99(mu, I, W2)
        return {'COUNTS': counts, 'I': I, 'MU': mu, 'W2': W2, 'MDP_99': mdp}

    def energy_binning(self):
        """Return the underlying energy binning in the form of a numpy array.

        Note the array has one more element than the underlying ENERG_LO and
        ENERG_HI arrays, i.e., if the cube is binned in two energy layers,
        say 2--4 keV and 4--8 keV, the energy binning is [2. 4. 8.].
        """
        return numpy.append(self.ENERG_LO, self.ENERG_HI[-1])

    def num_energy_layers(self):
        """Return the number of energy layers in the cube.
        """
        return len(self.ENERG_LO)

    def energy_range(self, energy_layer):
        """Return the energy range, i.e., (emin, emax) in keV, for a given
        energy layer.
        """
        return self.ENERG_LO[energy_layer], self.ENERG_HI[energy_layer]

    def energy_label(self, energy_layer):
        """Return a text label corresponding to the energy layer (e.g, to
        identify the energy range on a plot).
        """
        return '%.2f--%.2f keV' % self.energy_range(energy_layer)

    def _plot_layer(self, data, energy_layer, **kwargs):
        """Delegated function to make a plot of one of the underlying image
        extensions in *one* energy layer of the cube.

        This essentially calculates the proper slice of the input data and all
        the necessary arguments to be passed to xFITSImageBase for doing the
        actual plot.

        Note this function is *not* creating any matplotlib canvas---it is left
        to the caller to do that.

        Arguments
        ---------
        data : array
            The underlying data array. This can be the array corresponding to
            one of the underlying image extensions of the cube, or any
            array of the same shape.

        energy_layer : int
            The identifier of the energy layer in the cube.

        kwargs : dict
            Additional keyword arguments being passed to the underlying
            xFITSImageBase.make_plot() call.
        """
        assert energy_layer in range(0, self.num_energy_layers())
        data = data[energy_layer, :, :]
        slices = ('x', 'y', energy_layer)
        xFITSImageBase.make_plot(data, self.wcs, slices=slices, **kwargs)
        xFITSImageBase.add_label(self.energy_label(energy_layer), y=0.92)

    def _plot_base(self, data, energy_layers, mask, figname, prefix=None,
                   zlabel=None, post_plot_hook=None, **kwargs):
        """Delegated base function to plot multiple energy layers of the same
        image extension.
        """
        # pylint: disable=arguments-differ
        # Do some magic with the energy_layers argument.
        if energy_layers is None:
            energy_layers = list(range(0, self.num_energy_layers()))
        elif isinstance(energy_layers, int):
            energy_layers = [energy_layers]
        # If we are passing a mask, we create a copy of the data array (not to)
        # mess up with the underlying data, and explicitely set to zero all the
        # values that are not in the mask.
        if mask is not None:
            data = data.copy()
            data[numpy.logical_not(mask)] = 0.
        # Loop over the energy layers.
        if zlabel is None:
            zlabel = figname
        kwargs.setdefault('zlabel', zlabel)
        for energy_layer in energy_layers:
            figure_name = '%s %s' % (figname, self.energy_label(energy_layer))
            if prefix is not None:
                figure_name = '%s %s' % (prefix, figure_name)
            plt.figure(figure_name)
            self._plot_layer(data, energy_layer, **kwargs)
            if post_plot_hook is not None:
                post_plot_hook(energy_layer, mask)

    def plot_mdp_map(self, energy_layers=None, prefix=None, **kwargs):
        """Plot the MDP map.
        """
        self._plot_base(self.MDP_99, energy_layers, None, 'MDP', prefix, 'MDP (99% CL)', **kwargs)

    def plot(self, prefix=None):
        """Plot the data.
        """
        self.plot_mdp_map(prefix=prefix)

    def __str__(self):
        """String formatting.
        """
        text = xBinnedFileBase.__str__(self)
        for ext_name in self.EXTENSIONS_NAMES:
            val = ('%s' % self.__getattr__(ext_name)).replace('\n', '')
            text += '\n%s: %s' % (ext_name, val)
        text += '\nEnergy bounds:'
        for emin, emax in zip(self.ENERG_LO, self.ENERG_HI):
            text += '\n- %.3f--%.3f keV' % (emin, emax)
        return text



class xEventBinningMDPMAP(xEventBinningMDPMAPCUBE):

    """Class for MDPMAP binning.
    """

    INTENT = 'MDP map in sky coordinates'
    SUPPORTED_KWARGS = ['mc', 'acceptcorr', 'weights', 'weightcol', 'grayfilter',
        'emin', 'emax'] + xEventBinningBase._image_wcs_kwargs()

    def process_kwargs(self):
        """Overloaded method.

        Note we have to call the method of the base class closest in the
        inheritance hierarchy, to make sure that all the ingredients for the
        WCS object generation are properly setup.
        """
        xEventBinningMDPMAPCUBE.process_kwargs(self)
        self.set('ebins', 1)



class xEventBinningPMAPCUBE(xEventBinningMDPMAPCUBE):

    """Class for PMAPCUBE binning.
    """

    INTENT = 'polarization map cube in sky coordinates and energy layers'
    SUPPORTED_KWARGS = xEventBinningMDPMAPCUBE.SUPPORTED_KWARGS
    EXTENSIONS_NAMES = xBinTableHDUPCUBE.POL_COL_NAMES

    @staticmethod
    def _binned_data(stokes_analysis, *args, **kwargs):
        """Retrieve the necessary data from the underlying xStokesAnalysis object.
        """
        return stokes_analysis.polarization_map_cube(*args, **kwargs)



class xEventBinningPMAP(xEventBinningPMAPCUBE):

    """Class for PMAP binning.
    """

    INTENT = 'polarization map in sky coordinates'
    SUPPORTED_KWARGS = xEventBinningMDPMAP.SUPPORTED_KWARGS

    def process_kwargs(self):
        """Overloaded method.

        Note we have to call the method of the base class closest in the
        inheritance hierarchy, to make sure that all the ingredients for the
        WCS object generation are properly setup.
        """
        xEventBinningPMAPCUBE.process_kwargs(self)
        self.set('ebins', 1)



class xBinnedPolarizationMapCube(xBinnedMDPMapCube):

    """Read-mode interface to a PMAPCUBE FITS file.
    """

    EXTENSIONS_NAMES = xEventBinningPMAPCUBE.EXTENSIONS_NAMES

    def __recalculate(self):
        """
        """
        # Recalculate the normalized Stokes parameters, and propagate the errors.
        args = self.I, self.Q, self.U, self.MU, self.W2
        self.QN, self.UN, self.I_ERR, self.Q_ERR, self.U_ERR, self.QN_ERR, self.UN_ERR, \
            self.QUN_COV, self.P_VALUE, self.CONFID, self.SIGNIF = \
                xStokesAnalysis.calculate_stokes_errors(*args)
        # Recalculate the polarization degree and angle.
        self.PD, self.PD_ERR, self.PA, self.PA_ERR = \
            xStokesAnalysis.calculate_polarization(*args, degrees=True)
        self.MDP_99 = xStokesAnalysis.calculate_mdp99(self.MU, self.I, self.W2)

    def __iadd__(self, other):
        """Overloaded method for binned data addition.
        """
        self._check_iadd(other, ('Q', 'U'))
        xBinnedMDPMapCube.__iadd__(self, other)
        self.Q += other.Q
        self.U += other.U
        self.__recalculate()
        return self

    def convolve(self, kernel):
        """Convolve the polarization cube with a generic binned kernel.
        """
        for layer in range(self.num_energy_layers()):
            self.I[layer] = scipy.signal.convolve2d(self.I[layer], kernel, mode='same')
            self.Q[layer] = scipy.signal.convolve2d(self.Q[layer], kernel, mode='same')
            self.U[layer] = scipy.signal.convolve2d(self.U[layer], kernel, mode='same')
            self.W2[layer] = scipy.signal.convolve2d(self.W2[layer], kernel, mode='same')
            self.MU[layer] = scipy.signal.convolve2d(self.MU[layer], kernel, mode='same') /\
                kernel.sum()
        self.__recalculate()

    def align(self):
        """Align the polarization direction to the radial direction on a bin by
        bin basis.

        Warning
        -------
        This will need to be reviewed when we fix issue #597.
        """
        x0, y0, _ = self.wcs.wcs.crpix
        _, dec0, _ = self.wcs.wcs.crval
        shape = self.map_shape()
        x, y = numpy.indices(shape)
        dx = (x - x0) * numpy.cos(numpy.radians(dec0))
        dy = y - y0
        phi0 = numpy.arctan2(dy, -dx)
        q0 = xStokesAnalysis.stokes_q(phi0, weights=None)
        u0 = xStokesAnalysis.stokes_u(phi0, weights=None)
        self.Q, self.U = align_stokes_parameters(self.Q, self.U, q0, u0)

    def radial_profile(self, n=10, layer=0):
        """Create a radial polarization profile of the source.

        Warning
        -------
        This needs to be cleaned up and properly documented.
        """
        x0, y0, _ = self.wcs.wcs.crpix
        shape = self.map_shape()
        _, pixel_size, _ = self.wcs.wcs.cdelt
        x, y = numpy.indices(shape)
        profile = []
        for r0 in range(n, shape[0] // 2 - n):
            r = numpy.sqrt((x - x0)**2. + (y - y0)**2)
            mask = numpy.logical_and(r > r0 - n, r < r0 + n)
            data = self.sum_pixels(mask, layer)
            data['RADIUS'] = r0 * pixel_size * 60.
            profile.append(data)
        return profile

    def sum_pixels(self, spatial_mask, energy_layer=0):
        """Sum the relevant quantities over a given spatial mask and energy slice.
        """
        args = spatial_mask, energy_layer
        data = xBinnedMDPMapCube.sum_pixels(self, *args)
        data['Q'] = self._sum_image_pixels(self.Q, *args)
        data['U'] = self._sum_image_pixels(self.U, *args)
        args = data['I'], data['Q'], data['U'], data['MU'], data['W2']
        pd, pd_err, pa, pa_err = xStokesAnalysis.calculate_polarization(*args, degrees=True)
        data.update({'PD': pd, 'PD_ERR': pd_err, 'PA': pa, 'PA_ERR': pa_err})
        qn, un, i_err, q_err, u_err, qn_err, un_err, qun_cov, p_value, confid, signif = \
            xStokesAnalysis.calculate_stokes_errors(*args)
        data.update({'QN': qn, 'UN': un, 'I_ERR': i_err, 'Q_ERR': q_err, 'U_ERR': u_err,
            'QN_ERR': qn_err, 'UN_ERR': un_err, 'QUN_COV': qun_cov, 'P_VALUE': p_value,
            'CONFID': confid, 'SIGNIF': signif})
        return data

    def calculate_significance_mask(self, num_sigma=2., intensity_percentile=0.):
        """Utitlity function to calculate a mask on the underlying pixel cube
        where the polarization degree and angle can be reliably plotted.

        The calculation is based on two different metrics:

        * the value of I in each given pixel, compared with a given percentile
          of the overall I distribution across the cube;
        * the number of sigma that the measured polarization degree differs
          from either 0 or 1 (if the measurement is too close to 0 or 1 and
          misses the physical bounds by less than a few error bars, then it
          is effectively not a measurement).

        Note the logic is positive---i.e., we're passing the mask identifying
        the *good* pixels. This is typically negated downstream to set to
        zero all the other pixels.

        Arguments
        ---------
        num_sigma : float
            The number of standard deviations representing the significance of
            the measurement of the polarization degree.

        intensity_percentile : float
            the threshold on the I value for the pixels, based on the percentile
            of the I distribution across all the pixels.
        """
        i_mask = self.I >= numpy.percentile(self.I[self.I > 0], intensity_percentile)
        pd_mask_lo = self.PD >= num_sigma * self.PD_ERR
        pd_mask_hi = self.PD <= 1. - num_sigma * self.PD_ERR
        # We need this, in case we set the errors to zero for bins with no counts.
        err_mask = numpy.logical_and(self.PD_ERR > 0., self.PA_ERR > 0.)
        return numpy.logical_and.reduce((i_mask, pd_mask_lo, pd_mask_hi, err_mask))

    def plot_stokes_parameters(self, energy_layers=None, prefix=None, **kwargs):
        """Plot the Stokes parameters.
        """
        self._plot_base(self.I, energy_layers, None, 'Stokes I', prefix, **kwargs)
        self._plot_base(self.Q, energy_layers, None, 'Stokes Q', prefix, **kwargs)
        self._plot_base(self.U, energy_layers, None, 'Stokes U', prefix, **kwargs)

    def plot_significance(self, energy_layers=None, prefix=None, **kwargs):
        """Plot the significance.
        """
        self._plot_base(self.SIGNIF, energy_layers, None, 'Significance', prefix, **kwargs)

    def normalized_stokes_parameters(self):
        """Calculate the normalized Stokes parameters, i.e., Q/I and U/I.
        """
        q = xModelStokesParameters.normalize(self.Q, self.I)
        u = xModelStokesParameters.normalize(self.U, self.I)
        return q, u

    def plot_normalized_stokes_parameters(self, energy_layers=None, num_sigma=1.,
                                          intensity_percentile=0.05, prefix=None,
                                          **kwargs):
        """Plot the normalized Stokes parameters.
        """
        q, u = self.normalized_stokes_parameters()
        mask = self.calculate_significance_mask(num_sigma, intensity_percentile)
        self._plot_base(q, energy_layers, mask, 'Stokes Q/I', prefix, **kwargs)
        self._plot_base(u, energy_layers, mask, 'Stokes U/I', prefix, **kwargs)

    def _calculate_sky_grid(self):
        """Calculate the sky grid corresponding to the center of the pixels
        in the underlying map extensions.

        Note that this could be cached and calculated exactly once, but I
        doubt it would make any practical difference.
        """
        nx, ny = self.map_shape()
        x, y = numpy.meshgrid(numpy.arange(nx), numpy.arange(ny))
        ra, dec, _ = self.wcs.wcs_pix2world(x, y, [0], 0)
        return ra, dec

    def _overlay_arrows(self, energy_layer, mask=None, **kwargs):
        """Method to overlay the arrows of the polarization information.
        """
        x, y = self._calculate_sky_grid()
        args = self.PD[energy_layer], self.PA[energy_layer]
        dx, dy = xModelStokesParameters.pdpa_to_xy(*args, degrees=True)
        if mask is not None:
            mask = mask[energy_layer]
            x, y, dx, dy = x[mask], y[mask], dx[mask], dy[mask]
        plot_arrows((x, y), (dx, dy), **kwargs)

    def plot_polarization_degree(self, energy_layers=None, num_sigma=2.,
                                 arrows=True, prefix=None, **kwargs):
        """Plot the polarization degree.
        """
        mask = self.calculate_significance_mask(num_sigma, 0.)
        if arrows:
            hook = self._overlay_arrows
        else:
            hook = None
        self._plot_base(self.PD, energy_layers, mask, 'Polarization degree',
                        prefix, post_plot_hook=hook, **kwargs)

    def plot_polarization_angle(self, energy_layers=None, num_sigma=2.,
                                prefix=None, **kwargs):
        """Plot the polarization angle.
        """
        mask = self.calculate_significance_mask(num_sigma, 0.)
        self._plot_base(self.PA, energy_layers, mask, 'Polarization angle',
                        prefix, 'Polarization angle [$^\\circ$]', **kwargs)

    def save_to_ds9(self, output_file_path, num_sigma=2., intensity_percentile=0.,
        line_color='white'):
        """Save the polarization arrows to a ds9 region.
        """
        # pylint: disable=line-too-long
        _x, _y = self._calculate_sky_grid()
        sigma_mask = self.calculate_significance_mask(num_sigma, intensity_percentile)
        for energy_layer in range(self.num_energy_layers()):
            pd, pa = self.PD[energy_layer], self.PA[energy_layer]
            mask = sigma_mask[energy_layer]
            x = _x[mask]
            y = _y[mask]
            masked_pd = pd[mask]
            masked_pa = pa[mask]
            file_path = '%s_%d.reg' % (output_file_path, energy_layer)
            logger.info('Saving the polarization values ds9 region file to %s...', file_path)
            with open(file_path, 'w') as ds9_region_file:
                ds9_region_file.write('# Region file for ixpe polarization degree and angle at %s sigma significance\n' % num_sigma)
                ds9_region_file.write('global color=%s dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n fk5\n' % line_color)
                # Multiplying the pd value by 60 because the vector length in the
                # ds9 region is in arcseconds.
                v_length = arcmin_to_arcsec(masked_pd)
                for x_pos, y_pos, vl, pa in zip(x, y, v_length, masked_pa):
                    ds9_region_file.write('#vector(%s,%s,%s",%s) vector=0\n' %\
                        (x_pos, y_pos, vl, pa))

    def plot(self, prefix=None):
        """Plot everything.
        """
        self.plot_stokes_parameters(prefix=prefix)
        self.plot_normalized_stokes_parameters(prefix=prefix)
        self.plot_polarization_degree(prefix=prefix)
        self.plot_polarization_angle(prefix=prefix)

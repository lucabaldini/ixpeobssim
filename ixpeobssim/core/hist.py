#!/usr/bin/env python
#
# Copyright (C) 2018--2020, the ixpeobssim team.
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

"""Histogram facilities.
"""

from __future__ import print_function, division

import numbers

from astropy.io import fits
import matplotlib
import numpy
import scipy.stats
import scipy.signal

from ixpeobssim.core.fitting import fit_histogram, USE_ABSOLUTE_SIGMA
from ixpeobssim.instrument.gpd import gpd_map_binning, GPD_PHYSICAL_HALF_SIDE_X,\
    GPD_PHYSICAL_HALF_SIDE_Y
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, draggable_colorbar
from ixpeobssim.utils.logging_ import logger


# pylint: disable=invalid-name, too-many-arguments, attribute-defined-outside-init, too-many-instance-attributes


class xHistogramBase:

    """Base class for an n-dimensional histogram.

    This interface to histograms is profoundly different for the minimal
    numpy/matplotlib approach, where histogramming methods return bare
    vectors of bin edges and counts. The main underlying ideas are

    * we keep track of the bin contents, the bin entries and the sum of the
      weights squared (the latter for the purpose of calculating the errors);
    * we control the axis label and the plotting styles;
    * we provide two separate interfaces, fill() and set_content(), to
      fill the histogram from either unbinned or binned data;
    * we support the basic arithmetics (addition, subtraction and multiplication
      by a scalar);
    * we support full data persistence (I/O) in FITS format.

    Note that this base class is not meant to be instantiated directly, and
    the interfaces to concrete histograms of specific dimensionality are
    defined in the sub-classes.

    Parameters
    ----------
    binning : n-tuple of array
        the bin edges on the different axes.

    labels : n-tuple of strings
        the text labels for the different axes.
    """

    PLOT_OPTIONS = dict()

    def __init__(self, binning, labels):
        """Constructor.
        """
        assert len(labels) == len(binning) + 1
        self.binning = tuple(binning)
        self.labels = list(labels)
        self.shape = tuple(len(bins) - 1 for bins in self.binning)
        self.num_axes = len(self.binning)
        self.content = numpy.zeros(shape=self.shape, dtype=float)
        self.entries = numpy.zeros(shape=self.shape, dtype=float)
        self.sumw2 = numpy.zeros(shape=self.shape, dtype=float)

    def fill(self, *data, weights=None):
        """Fill the histogram from unbinned data.

        Note this method is returning the histogram instance, so that the function
        call can be chained.
        """
        if weights is None:
            weights = numpy.ones(data[0].shape, dtype=float)
        elif isinstance(weights, numbers.Number):
            weights = numpy.full(data[0].shape, weights, dtype=float)
        data = numpy.vstack(data).T
        content, _ = numpy.histogramdd(data, bins=self.binning, weights=weights)
        entries, _ = numpy.histogramdd(data, bins=self.binning)
        sumw2, _ = numpy.histogramdd(data, bins=self.binning, weights=weights**2.)
        self.content += content
        self.entries += entries
        self.sumw2 += sumw2
        return self

    def set_content(self, content, entries=None, errors=None):
        """Set the bin contents programmatically from binned data.

        Note this method is returning the histogram instance, so that the function
        call can be chained.
        """
        assert content.shape == self.shape
        self.content = content
        if entries is not None:
            assert entries.shape == self.shape
            self.entries = entries
        if errors is not None:
            self.set_errors(errors)
        return self

    def errors(self):
        """Return the bin errors.
        """
        return numpy.sqrt(self.sumw2)

    def set_errors(self, errors):
        """
        """
        assert errors.shape == self.shape
        self.sumw2 = errors**2.

    def bin_centers(self, axis=0):
        """Return the bin centers for a specific axis.
        """
        return 0.5 * (self.binning[axis][1:] + self.binning[axis][:-1])

    def bin_widths(self, axis=0):
        """Return the bin widths for a specific axis.
        """
        return numpy.diff(self.binning[axis])

    def _axis_complement(self, axis=0):
        """Return a tuple of integers identifying all the histogram axes
        but the one passed as an arguments. If, e.g, num_axes is 3,
        calling _axis_complement(1) will return (0, 2).

        This is used to sum and/or average on all the histogram dimensions
        except for a given one, see, e.g., the mean() and rms() methods below.
        """
        return tuple(ax for ax in range(self.num_axes) if ax != axis)

    def mean(self, axis=0):
        """Calculate the (binned) mean eastimate along a given axis.
        """
        content = numpy.sum(self.content, axis=self._axis_complement(axis))
        return numpy.sum(self.bin_centers(axis) * content) / self.sum()

    def rms(self, axis=0):
        """Calculate the (binned) rms eastimate along a given axis.

        .. warning::
           Mind this is not using anything clever (e.g. the Welford algorithm)
           and is not particularly numerically stable.
        """
        bin_centers = self.bin_centers(axis)
        content = numpy.sum(self.content, axis=self._axis_complement(axis))
        _sum = numpy.sum(bin_centers * content) / self.sum()
        _sum2 = numpy.sum(bin_centers**2. * content) / self.sum()
        return numpy.sqrt(_sum2 - _sum**2.)

    @staticmethod
    def bisect(binning, values, side='left'):
        """Return the indices corresponding to a given array of values for a
        given binning.
        """
        return numpy.searchsorted(binning, values, side) - 1

    def find_bin(self, *coords):
        """Find the bin corresponding to a given set of "physical" coordinates
        on the histogram axes.

        This returns a tuple of integer indices that can be used to address
        the histogram content.
        """
        return tuple([self.bisect(binning, value) for binning, value in zip(self.binning, coords)])

    def find_bin_value(self, *coords):
        """Find the histogram content corresponding to a given set of "physical"
        coordinates on the histogram axes.
        """
        return self.content[self.find_bin(*coords)]

    def num_entries(self):
        """Return the number of entries in the histogram.
        """
        return self.entries.sum()

    def sum(self, axis=None):
        """return the sum of weights in the histogram.
        """
        return self.content.sum(axis)

    def empty_copy(self):
        """Create an empty copy of a histogram.
        """
        return self.__class__(*self.binning, *self.labels)

    def copy(self):
        """Create a full copy of a histogram.
        """
        hist = self.empty_copy()
        hist.set_content(self.content.copy(), self.entries.copy())
        return hist

    def __add__(self, other):
        """Histogram addition.
        """
        hist = self.empty_copy()
        hist.set_content(self.content + other.content,
                         self.entries + other.entries,
                         numpy.sqrt(self.sumw2 + other.sumw2))
        return hist

    def __sub__(self, other):
        """Histogram subtraction.
        """
        hist = self.empty_copy()
        hist.set_content(self.content - other.content,
                         self.entries + other.entries,
                         numpy.sqrt(self.sumw2 + other.sumw2))
        return hist

    def __mul__(self, value):
        """Histogram multiplication by a scalar.

        Args
        ----
        value : array_like
            The scale factor for the multiplication---must be either a scalar
            or an array of the same shape of the histogram content.
        """
        hist = self.empty_copy()
        hist.set_content(self.content * value, self.entries, self.errors() * value)
        return hist

    def __rmul__(self, value):
        """Histogram multiplication by a scalar.
        """
        return self.__mul__(value)

    def set_axis_label(self, axis, label):
        """Set the label for a given axis.
        """
        self.labels[axis] = label

    def _plot(self, **kwargs):
        """No-op plot() method, to be overloaded by derived classes.
        """
        raise NotImplementedError('_plot() not implemented for %s' % self.__class__.__name__)

    def plot(self, **kwargs):
        """Plot the histogram.
        """
        for key, value in self.PLOT_OPTIONS.items():
            kwargs.setdefault(key, value)
        self._plot(**kwargs)
        setup_gca(xmin=self.binning[0][0], xmax=self.binning[0][-1],
                  xlabel=self.labels[0], ylabel=self.labels[1])

    @staticmethod
    def label_keyword(axis):
        """Header keyword for the axis labels.
        """
        return 'LABEL%d' % axis

    @staticmethod
    def entries_hdu_name():
        """Extension name for the bin entries.
        """
        return 'ENTRIES'

    @staticmethod
    def sumw2_hdu_name():
        """Extension name for the bin entries.
        """
        return 'SUMW2'

    @staticmethod
    def binning_hdu_name(axis):
        """Extension name for the axis binnings.
        """
        return 'BINNING%d' % axis

    @staticmethod
    def binning_col_name():
        """Column name for the axis binnings.
        """
        return 'EDGES'

    @classmethod
    def from_file(cls, file_path):
        """ Load the histogram from a FITS file.

        Note that we transpose the image data back at read time---see the
        comment in the save() method.
        """
        logger.info('Loading histogram from %s...', file_path)
        with fits.open(file_path) as hdu_list:
            content = hdu_list['Primary'].data.T
            entries = hdu_list[cls.entries_hdu_name()].data.T
            errors = numpy.sqrt(hdu_list[cls.sumw2_hdu_name()].data.T)
            num_axes = len(content.shape)
            edges = [hdu_list[cls.binning_hdu_name(i)].data[cls.binning_col_name()] \
                for i in range(num_axes)]
            labels = [hdu_list['Primary'].header[cls.label_keyword(i)] for i in range(num_axes + 1)]
        hist = cls(*edges, *labels)
        hist.set_content(content, entries, errors)
        return hist

    def save(self, file_path, overwrite=True, **header_keywords):
        """ Save the histogram (with edges) to a FITS file.

        Note that all the image data are transposed so that the thing can be
        correctly visualized with standard FITS viewers---at least in two
        dimensions.
        """
        logger.info('Saving histogram to %s...', file_path)
        # Save the histogram content in a series of images.
        hdu_list = [
            fits.PrimaryHDU(self.content.T),
            fits.ImageHDU(self.entries.T, name=self.entries_hdu_name()),
            fits.ImageHDU(self.sumw2.T, name=self.sumw2_hdu_name())
        ]
        for key, value in header_keywords.items():
            hdu_list[0].header.set(key, value)
        for i, label in enumerate(self.labels):
            hdu_list[0].header.set(self.label_keyword(i), label)
        for i, binning in enumerate(self.binning):
            col = fits.Column(name=self.binning_col_name(), array=binning, format='E')
            hdu = fits.BinTableHDU.from_columns([col])
            hdu.name = self.binning_hdu_name(i)
            hdu_list.append(hdu)
        hdu_list = fits.HDUList(hdu_list)
        hdu_list.writeto(file_path, overwrite=overwrite)



class xHistogram1d(xHistogramBase):

    """Container class for one-dimensional histograms.
    """

    PLOT_OPTIONS = dict(lw=1.25, alpha=0.4, histtype='stepfilled')

    def __init__(self, xbins, xlabel='', ylabel='Entries/bin'):
        """Constructor.
        """
        xHistogramBase.__init__(self, (xbins, ), [xlabel, ylabel])

    def _plot(self, **kwargs):
        """Overloaded make_plot() method.
        """
        plt.hist(self.bin_centers(0), self.binning[0], weights=self.content, **kwargs)

    def errorbar_data(self):
        """Return the x, y, dy arrays that can be used to build a scatter plot
        (with errors) from the histogram.
        """
        return self.bin_centers(0), self.content, self.errors()

    def errorbar(self, **kwargs):
        """Plot the histogram as a scatter plot.
        """
        kwargs.setdefault('fmt', 'o')
        x, y, dy = self.errorbar_data()
        plt.errorbar(x, y, dy, **kwargs)
        setup_gca(xlabel=self.labels[0], ylabel=self.labels[1])

    def scatter_plot(self):
        """Turn the histogram into a scatter plot.

        .. warning::
           This is to be removed.
        """
        return xScatterPlot(self.bin_centers(0), self.content, xlabel=self.labels[0])

    def gaussian_kde_smooth(self, bandwidth=2.):
        """Create a copy of the histogram where the weights are smoothed with a
        gaussian kernel density estimatore.

        Args
        ----
        bandwidth : float
            The sigma of the gaussian kernel, in (fractional) number of bins.
        """
        n = int(5. * bandwidth + 0.5)
        kernel_domain = numpy.linspace(-n, n, 2 * n + 1)
        kernel = scipy.stats.norm(0, bandwidth).pdf(kernel_domain)
        hist = self.copy()
        hist.content = numpy.convolve(self.content, kernel, 'same')
        return hist

    def fit(self, model, p0=None, sigma=None, xmin=-numpy.inf, xmax=numpy.inf,
            absolute_sigma=USE_ABSOLUTE_SIGMA, check_finite=True, method=None,
            verbose=True, **kwargs):
        """Fit the histogram with a model.
        """
        return fit_histogram(model, self, p0, sigma, xmin, xmax, absolute_sigma,
                             check_finite, method, verbose, **kwargs)



class xScatterPlot:

    """Small class encapsulating a scatter plot.

    Technically speaking, this would not belong here, as a scatter plot is not
    strictly related to any of the histogram classes, but a 1-dimensional
    histogram with errors can techincally be turned into a scatter plot,
    so we introduce the concept here for completeness.

    .. warning::
        Consider removing this class.
    """

    def __init__(self, x, y, dy=None, dx=None, xlabel=None, ylabel=None):
        """Constructor.
        """
        self.x = x
        self.y = y
        self.dy = dy
        self.dx = dx
        self.labels = (xlabel, ylabel)

    def plot(self, **kwargs):
        """Plot the scatter plot.
        """
        kwargs.setdefault('fmt', 'o')
        kwargs.setdefault('markeredgecolor', 'black')
        plt.errorbar(self.x, self.y, self.dy, self.dx, **kwargs)
        setup_gca(xlabel=self.labels[0], ylabel=self.labels[1])



class xHistogram2d(xHistogramBase):

    """Container class for two-dimensional histograms.
    """

    PLOT_OPTIONS = dict(cmap=plt.get_cmap('hot'))

    def __init__(self, xbins, ybins, xlabel='', ylabel='', zlabel='Entries/bin'):
        """Constructor.
        """
        xHistogramBase.__init__(self, (xbins, ybins), [xlabel, ylabel, zlabel])

    def _plot(self, logz=False, **kwargs):
        """Overloaded make_plot() method.
        """
        x, y = (v.flatten() for v in numpy.meshgrid(self.bin_centers(0), self.bin_centers(1)))
        bins = self.binning
        w = self.content.T.flatten()
        if logz:
            # Hack for a deprecated functionality in matplotlib 3.3.0
            # Parameters norm and vmin/vmax should not be used simultaneously
            # If logz is requested, we intercent the bounds when created the norm
            # and refrain from passing vmin/vmax downstream.
            vmin = kwargs.pop('vmin', None)
            vmax = kwargs.pop('vmax', None)
            kwargs.setdefault('norm', matplotlib.colors.LogNorm(vmin, vmax))
        plt.hist2d(x, y, bins, weights=w, **kwargs)
        colorbar = draggable_colorbar()
        if self.labels[2] is not None:
            colorbar.set_label(self.labels[2])

    def gaussian_kde_smooth(self, bandwidth=(2., 2.)):
        """Create a copy of the histogram where the weights are smoothed with a
        gaussian kernel density estimatore.

        .. warning::
           This is essentially untested.

        Args
        ----
        bandwidth : 2-element tuple float
            The sigma of the gaussian kernel, in (fractional) number of bins.
        """
        nx, ny = [int(5. * bw + 0.5) for bw in bandwidth]
        dx = numpy.linspace(-nx, nx, 2 * nx + 1)
        dy = numpy.linspace(-ny, ny, 2 * ny + 1)
        kernel_domain = numpy.dstack(numpy.meshgrid(dx, dy))
        sigma = numpy.array([[bandwidth[0], 0.], [0., bandwidth[1]]])
        kernel = scipy.stats.multivariate_normal((0., 0.), sigma).pdf(kernel_domain)
        hist = self.copy()
        hist.content = scipy.signal.convolve2d(self.content, kernel, 'same')
        return hist

    def hslice(self, bin_):
        """Return the horizontal slice for a given bin.
        """
        hist = xHistogram1d(self.binning[0], self.labels[0])
        hist.set_content(self.content[:, bin_], self.entries[:, bin_])
        return hist

    def hslices(self):
        """Return a list of all the horizontal slices.
        """
        return tuple(self.hslice(bin_) for bin_ in range(self.shape[1]))

    def hbisect(self, y):
        """Return the horizontal slice corresponding to a given y value.
        """
        return self.hslice(self.bisect(self.binning[1], y))

    def vslice(self, bin_):
        """Return the vertical slice for a given bin.
        """
        hist = xHistogram1d(self.binning[1], self.labels[1])
        hist.set_content(self.content[bin_], self.entries[bin_])
        return hist

    def vslices(self):
        """Return a list of all the vertical slices.
        """
        return tuple(self.vslice(bin) for bin in range(self.shape[0]))

    def vbisect(self, x):
        """Return the vertical slice corresponding to a given y value.
        """
        return self.vslice(self.bisect(self.binning[0], x))



class xModulationCube2d(xHistogram2d):

    """Specialized class for modulations cubes.
    """

    def __init__(self, xbins, ybins):
        """Constructor.
        """
        labels = ['Energy [keV]', '$\\phi$ [rad]', 'Entries/bin']
        xHistogram2d.__init__(self, xbins, ybins, *labels)



class xGpdMap2d(xHistogram2d):

    """2-dimensional GPD map.
    """

    def __init__(self, nside=10, zlabel='Entries/bin'):
        """Constructor.
        """
        edges = gpd_map_binning(GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_Y, nside)
        labels = ['x [mm]', 'y [mm]', zlabel]
        xHistogram2d.__init__(self, *edges, *labels)



class xHistogram3d(xHistogramBase):

    """Container class for three-dimensional histograms.
    """

    def __init__(self, xbins, ybins, zbins, xlabel='', ylabel='', zlabel='',
                 wlabel='Entries/bin'):
        """Constructor.
        """
        xHistogramBase.__init__(self, (xbins, ybins, zbins), [xlabel, ylabel, zlabel, wlabel])



class xGpdMap3d(xHistogram3d):

    """Three-dimensional histogram where the first two axes represent the
    GPD active area in Physical coordinates.
    """

    def __init__(self, nside, zbins, zlabel='', wlabel='Entries/bin'):
        """Constructor.
        """
        edges = gpd_map_binning(GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_Y, nside)
        xHistogram3d.__init__(self, *edges, zbins, 'x [mm]', 'y [mm]', zlabel, wlabel)

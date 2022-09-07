#!/usr/bin/env python
#
# Copyright (C) 2020, the ixpeobssim team.
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

from __future__ import print_function, division

import numpy
from astropy.io import fits

from ixpeobssim.core.fitsio import xPrimaryHDU, xBinTableHDUBase
from ixpeobssim.core.fitsio import FITS_TO_NUMPY_TYPE_DICT
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.core.modeling import xGaussian
from ixpeobssim.core.fitting import fit_histogram
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, xStatBox
from ixpeobssim.utils.math_ import format_value_error



class xBenchmarkTable:

    """Simple class representing a flat table.

    This is meant to be a naive interface to save a flat table of values to a
    FITS file, with the additional convenience of being able to add the values
    by row, rather than preparing the columns externally.

    Note that, for each column, we add a fellow one meant to hold the
    corresponding error. Encapsulating this very logic into the class makes it
    easy to calculate quantities, such as the pulls, where the values and
    errors go hand in hand.

    Warning
    -------
    This is dumb, and is meant to. It's not nice-looking, nor optimized, and it
    does not even give you control on the cell format.
    """

    FITS_TYPE = 'E'
    DTYPE = FITS_TO_NUMPY_TYPE_DICT[FITS_TYPE]

    def __init__(self, col_names, creator='N/A', keywords=None, comments=None):
        """Constructor.
        """
        self.col_names = []
        for col_name in col_names:
            self.col_names += [col_name, self.error_label(col_name)]
        self.creator = creator
        self.keywords = keywords
        self.comments = comments
        self._specs = [(col_name, self.FITS_TYPE) for col_name in self.col_names]
        self._data = numpy.empty(shape=(0, len(self._specs)), dtype=self.DTYPE)

    @staticmethod
    def error_label(col_name):
        """Build the name for the error column associated with a generic value.
        """
        return '%s_ERR' % col_name

    def add_row(self, values):
        """Add a row to the table.
        """
        values = numpy.array(values, dtype=self.DTYPE)
        self._data = numpy.vstack((self._data, values))

    def writeto(self, file_path, overwrite=True):
        """Save the table to a FITS file.

        This creates subclasses xBinTableHDUBase under the hood and creates
        an instance of the new class, which is immediately written to file.
        """

        class __Table(xBinTableHDUBase):

            NAME = 'BENCHMARK'
            DATA_SPECS = self._specs

        cols = [self._data[:, i] for i in range(len(self.col_names))]
        _table = __Table(cols)
        _header = xPrimaryHDU(creator=self.creator, keywords=self.keywords,
                              comments=self.comments)
        _hdu_list = fits.HDUList([_header, _table])
        _hdu_list.info()
        logger.info('Writing benchmark table to %s...', file_path)
        _hdu_list.writeto(file_path, overwrite=overwrite)
        logger.info('Done.')



def load_benchmark_data(file_path):
    """Read back a FITS file.
    """
    logger.info('Loading benchmark data from %s...', file_path)
    with fits.open(file_path) as hdu_list:
        hdu_list.info()
        header = hdu_list[0].header
        data = hdu_list[1].data
    return header, data


def plot_pull_histogram(data, col_name, target=None, span_sigma=7.5, num_bins=100):
    """General-purpose method for plotting the pulls of a parameter estimate.
    """
    # Retrieve the values, the errors and the pulls.
    values = data[col_name]
    errors = data[xBenchmarkTable.error_label(col_name)]
    if target is None:
        target = values.mean()
    pulls = (values - target) / errors
    # Create the histogram.
    try:
        center = int(pulls.mean())
    except:
        center = 0.
    binning = numpy.linspace(center - span_sigma, center + span_sigma, num_bins)
    hist = xHistogram1d(binning, xlabel='Pull [$\\sigma$]').fill(pulls)
    hist.plot()
    # Fit the histogram.
    try:
        model = fit_histogram(xGaussian(), hist)
    except RuntimeError:
        return
    model.set_plotting_range(binning.min(), binning.max())
    model.plot()
    model.stat_box()
    # Add a custom stat box.
    mean = values.mean()
    mean_err = values.std(ddof=1) / numpy.sqrt(len(values))
    bias = (mean - target) / target
    box = xStatBox()
    box.add_entry(col_name)
    box.add_entry('Average', mean, mean_err)
    box.add_entry('Target', target)
    box.add_entry('Bias: %.2f%%' % (100. * bias))
    box.plot()
    setup_gca(grids=True)
    return hist, model


def plot_metrics_comparison(data, col_names, binning, labels=None, **kwargs):
    """Overlay a series of histograms for specific columns.
    """
    if labels is None:
        labels = col_names
    for col_name, label in zip(col_names, labels):
        col = data[col_name]
        hist = xHistogram1d(binning).fill(col)
        hist.plot(label=label)
    setup_gca(**kwargs)

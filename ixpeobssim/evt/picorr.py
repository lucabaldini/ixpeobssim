# Copyright (C) 2022, the ixpeobssim team.
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

"""Facilities for the PI correction (aka energy scaling).
"""

from __future__ import print_function, division

from astropy.io import fits
import numpy

from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, metplot, setup_gca
from ixpeobssim.utils.os_ import check_input_file
from ixpeobssim.utils.time_ import met_to_string


class xTrapezoidPdf:

    """Convenience class representing a trapeziod pdf for the purpose of
    randomizing the PI before the correction.

    This was setup to implement a second-order correction to the PI randomization,
    but it's not actually used.
    """

    @staticmethod
    def pdf(x, slope):
        """Probability density function.
        """
        return 1. + slope * x

    @staticmethod
    def cdf(x, slope):
        """Probability density function.
        """
        return 0.5 * slope * x**2. + x + 0.5 * (1 - 0.25 * slope)

    @staticmethod
    def ppf(q, slope):
        """Probability density function.
        """
        delta = 1. - 2. * slope * ((4. - slope) / 8. - q)
        return (numpy.sqrt(delta) - 1) / slope

    @staticmethod
    def rvs(slope, size=1):
        """
        """
        q = numpy.random.random(size)
        return xTrapezoidPdf.ppf(q, slope)

    @staticmethod
    def randomize(pi):
        """Randomize a vector of PI values.

        Args
        ----
        pi : array_like
            The array of PI values.
        """
        oversample  = 4
        binning = numpy.linspace(-0.5, 374.5, 95)
        hist = xHistogram1d(binning).fill(pi)
        hist.content /= oversample
        x = hist.bin_centers()
        y = hist.content
        spline = xInterpolatedUnivariateSpline(x, y, k=1)
        #spline.plot()
        m = spline(pi, 1) / spline(pi)
        print(m, m.max(), m.min())
        return pi + xTrapezoidPdf.rvs(m, pi.shape)



class xPulseInvariantCorrection:

    """Class describing a time-dependent linear correction for the
    pulse invariant, driven by a suitable FITS file.
    """

    EXT_NAME = 'PI_CALIBRATION'
    COL_NAMES = ('TIME', 'SLOPE', 'OFFSET')

    def __init__(self, file_path, k=1):
        """Constructor.
        """
        self._met, self._slope, self._offset = self._read_file(file_path)
        mask = numpy.logical_not(numpy.isnan(self._slope))
        self._slope_spline = xInterpolatedUnivariateSpline(self._met[mask], self._slope[mask], k=k)
        mask = numpy.logical_not(numpy.isnan(self._offset))
        self._offset_spline = xInterpolatedUnivariateSpline(self._met[mask], self._offset[mask], k=k)
        self.min_met = self._met.min()
        self.max_met = self._met.max()

    @staticmethod
    def _read_file(file_path):
        """Read the the correction data points from a file.
        """
        check_input_file(file_path, 'fits')
        logger.info('Reading PI correction coeffcients from %s...', file_path)
        col_names = xPulseInvariantCorrection.COL_NAMES
        with fits.open(file_path) as hdu_list:
            data = hdu_list[xPulseInvariantCorrection.EXT_NAME].data
            logger.info('%d row(s) found...', len(data))
            for col_name in col_names:
                num_nans = numpy.isnan(data[col_name]).sum()
                if num_nans > 0:
                    logger.warning('%s column has %d nan values.', col_name, num_nans)
            met, slope, offset = [data[col_name] for col_name in col_names]
        min_met = met.min()
        max_met = met.max()
        logger.info('Minimum MET = %.3f s (%s)', min_met, met_to_string(min_met))
        logger.info('Maximum MET = %.3f s (%s)', max_met, met_to_string(max_met))
        return met, slope, offset

    def __call__(self, met):
        """Return the interpolated slope and intercept at a given time or array of times.
        """
        num_left = (met < self.min_met).sum()
        if num_left > 0:
            logger.warning('%d MET values lower than the minimum correction MET...', num_left)
        num_right = (met > self.max_met).sum()
        if num_right > 0:
            logger.warning('%d MET values higher than the minimum correction MET...', num_right)
        if num_left or num_right:
            logger.warning('The PI correction splines will be extrapolated.')
        return self._slope_spline(met), self._offset_spline(met)

    def plot(self):
        """Plot the relevant quantities.
        """
        kwargs = dict(marker='o', lw=0)
        plt.figure('PI correction slope')
        metplot(self._met, self._slope, **kwargs)
        setup_gca(ylabel='Slope', grids=True)
        plt.figure('PI correction offset')
        metplot(self._met, self._offset, **kwargs)
        setup_gca(ylabel='Offset', grids=True)

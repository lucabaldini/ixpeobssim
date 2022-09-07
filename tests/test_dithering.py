#!/usr/bin/env python
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

from __future__ import print_function, division

import unittest
import sys

import numpy
import matplotlib.animation

from ixpeobssim.instrument.sc import dithering_pattern, spiral_dithering_pattern
from ixpeobssim.instrument.sc import triangular_wave, pow_triangular_wave
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.core.hist import xHistogram1d, xHistogram2d
from ixpeobssim.irf import load_psf, DEFAULT_IRF_NAME
from ixpeobssim.utils.units_ import degrees_to_arcmin

if sys.flags.interactive:
    plt.ion()



class TestDithering(unittest.TestCase):

    """Unit test for the dithering code.
    """

    @staticmethod
    def _figure(title):
        """
        """
        fig = plt.figure(title)
        plt.axis([-2., 2., -2., 2.])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('x [arcmin]')
        plt.ylabel('y [arcmin]')
        return fig

    @staticmethod
    def _hist2d(x, y, num_bins=50):
        """
        """
        binning = numpy.linspace(-2., 2., num_bins)
        return xHistogram2d(binning, binning).fill(x, y)

    def test_plain(self, num_bins=50):
        """
        """
        dithering = dithering_pattern()
        t = numpy.arange(0., 10000., 0.1)
        x, y = dithering(t)
        h2 = self._hist2d(x, y, num_bins)
        self._figure('Raw dithering')
        h2.plot()

    def test_psf_convolution(self, num_bins=100):
        """
        """
        amplitude = 1.6
        pa = 907.
        px = 101.
        py = 449
        dithering = dithering_pattern(amplitude, pa, px, py)
        t = numpy.arange(0., 40000., 0.002)
        x, y = dithering(t)
        psf = load_psf(DEFAULT_IRF_NAME, du_id=1)
        dx, dy = psf.delta(len(t))
        x += degrees_to_arcmin(dx)
        y += degrees_to_arcmin(dy)
        h2 = self._hist2d(x, y, num_bins)
        x, y = numpy.meshgrid(h2.bin_centers(0), h2.bin_centers(1))
        mask = (x**2. + y**2.) <= amplitude**2
        counts = h2.content[mask]
        binning = numpy.linspace(0., 2. * counts.max(), 100)
        h1 = xHistogram1d(binning).fill(counts.flatten())
        self._figure('PSF-convolved dithering')
        h2.plot()
        self._figure('PSF-convolved dithering with mask')
        h2.content[~mask] = 0.
        h2.plot()
        plt.figure('PSF-convolved dithering 1d')
        h1.plot()

    def test_psf_convolution_spiral(self, num_bins=100):
        """
        """
        amplitude = 1.6
        dithering = spiral_dithering_pattern(amplitude)
        t = numpy.arange(0., 40000., 0.002)
        x, y = dithering(t)
        psf = load_psf(DEFAULT_IRF_NAME, du_id=1)
        dx, dy = psf.delta(len(t))
        x += degrees_to_arcmin(dx)
        y += degrees_to_arcmin(dy)
        h2 = self._hist2d(x, y, num_bins)
        x, y = numpy.meshgrid(h2.bin_centers(0), h2.bin_centers(1))
        mask = (x**2. + y**2.) <= amplitude**2
        counts = h2.content[mask]
        binning = numpy.linspace(0., 2. * counts.max(), 100)
        h1 = xHistogram1d(binning).fill(counts.flatten())
        self._figure('PSF-convolved spiral dithering')
        h2.plot()
        self._figure('PSF-convolved spiral dithering with mask')
        h2.content[~mask] = 0.
        h2.plot()
        plt.figure('PSF-convolved spiral dithering 1d')
        h1.plot()

    def test_path(self):
        """
        """
        dithering = dithering_pattern()
        t = numpy.arange(0., 10001., 1)
        x, y = dithering(t)
        self._figure('Dithering path')
        plt.plot(x, y)
        plt.text(1.0, 1.5, '{} s'.format(t[-1]), fontsize=12)

    def test_path_spiral(self):
        """
        """
        dithering = spiral_dithering_pattern()
        t = numpy.arange(0., 5001., 1)
        x, y = dithering(t)
        self._figure('Spiral dithering path')
        plt.plot(x, y)
        plt.text(1.0, 1.5, '{} s'.format(t[-1]), fontsize=12)

    def test_triangular_wave(self, amplitude=1.6, period=1000.):
        """
        """
        plt.figure('Triangular wave')
        x = numpy.linspace(0., 5000., 1000)
        y = triangular_wave(x, amplitude, period)
        plt.plot(x, y, label='Simple triangular wave')
        y = pow_triangular_wave(x, amplitude, period, 0.5)
        plt.plot(x, y, label='Square root')
        plt.legend()


def test_animation(tstart=0., tstop=5000., tstep=2., interval=4):
    """Create a dithering animation.
    """
    # Prepare the data.
    dithering = dithering_pattern()
    t = numpy.arange(tstart, tstop, tstep)
    # Prepare the figure.
    TestDithering._figure('Dithering animation')
    # Prepare the artists.
    time_label = plt.text(1.0, 1.5, '', fontsize=12)
    line, = plt.gca().plot([], [], color='blue')

    def animate(step):
        """Small nested function doing the actual animation.
        """
        time_label.set_text('{} s'.format(t[step]))
        line.set_data(*dithering(t[:step]))
        return line, time_label

    # Shoot!
    args = plt.gcf(), animate, len(t)
    kwargs = dict(interval=interval, blit=True, repeat=False, save_count=50)
    return matplotlib.animation.FuncAnimation(*args, **kwargs)



if __name__ == '__main__':
    if False:
        animation = test_animation()
        #animation.save('dithering.gif', writer='imagemagick', fps=30)
        plt.show()
    else:
        unittest.main(exit=not sys.flags.interactive)

#!/usr/bin/env python
#
# Copyright (C) 2018--2021, the ixpeobssim team.
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
import os
import sys

import numpy

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.core.hist import xHistogram1d, xHistogram2d, xScatterPlot,\
    xModulationCube2d, xGpdMap2d
from ixpeobssim.utils.matplotlib_ import plt

numpy.random.seed(666)

if sys.flags.interactive:
    plt.ion()


class TestHist(unittest.TestCase):

    @classmethod
    def setUpClass(cls, sample_size=100000):
        """Set up the test.
        """
        cls.sample_size = sample_size
        cls.xedges = numpy.linspace(0, 1, 101)
        cls.yedges = numpy.linspace(2, 4, 101)
        cls.xmean = 0.5
        cls.xsigma = 0.1
        cls.x = numpy.random.normal(cls.xmean, cls.xsigma, size=cls.sample_size)
        cls.ymin = 2.
        cls.ymax = 4.
        cls.y = numpy.random.uniform(cls.ymin, cls.ymax, size=cls.sample_size)
        cls.h1 = xHistogram1d(cls.xedges, 'x [a. u.]')
        cls.h1.fill(cls.x)
        cls.h2 = xHistogram2d(cls.xedges, cls.yedges, 'x [a. u.]', 'y [a. u.]')
        cls.h2.fill(cls.x, cls.y)

    def test_plot(self):
        """
        """
        plt.figure('1d histogram')
        self.h1.plot()
        plt.figure('1d histogram scatter')
        self.h1.errorbar()
        plt.figure('2d histogram')
        self.h2.plot()

    def test_entries_unweighted(self):
        """
        """
        self.assertTrue(numpy.allclose(self.h1.entries, self.h1.content))
        self.assertTrue(numpy.allclose(self.h2.entries, self.h2.content))

    def test_weigths(self):
        """
        """
        h1 = self.h1.empty_copy()
        h1.fill(self.x, weights=numpy.full(self.x.shape, 2.))
        self.assertTrue(numpy.allclose(h1.entries, self.h1.entries))
        self.assertTrue(numpy.allclose(h1.content, 2 * self.h1.content))
        h2 = self.h2.empty_copy()
        h2.fill(self.x, self.y, weights=numpy.full(self.x.shape, 2.))
        self.assertTrue(numpy.allclose(h2.entries, self.h2.entries))
        self.assertTrue(numpy.allclose(h2.content, 2 * self.h2.content))

    def test_empty_copy(self):
        """Test the empty copy.
        """
        h1 = self.h1.empty_copy()
        self.assertEqual(h1.shape, self.h1.shape)
        for axis in range(h1.num_axes):
            self.assertTrue(numpy.allclose(h1.binning[axis], self.h1.binning[axis]))
        self.assertEqual(h1.labels, self.h1.labels)
        self.assertTrue(numpy.allclose(h1.entries, numpy.zeros(h1.shape, dtype=int)))
        self.assertTrue(numpy.allclose(h1.content, numpy.zeros(h1.shape, dtype=float)))
        h2 = self.h2.empty_copy()
        self.assertEqual(h2.shape, self.h2.shape)
        for axis in range(h2.num_axes):
            self.assertTrue(numpy.allclose(h2.binning[axis], self.h2.binning[axis]))
        self.assertEqual(h2.labels, self.h2.labels)
        self.assertTrue(numpy.allclose(h2.entries, numpy.zeros(h2.shape, dtype=int)))
        self.assertTrue(numpy.allclose(h2.content, numpy.zeros(h2.shape, dtype=float)))

    def test_copy(self):
        """Test copy.
        """
        h1 = self.h1.copy()
        self.assertTrue(numpy.allclose(h1.entries, self.h1.entries))
        self.assertTrue(numpy.allclose(h1.content, self.h1.content))
        h2 = self.h2.copy()
        self.assertTrue(numpy.allclose(h2.entries, self.h2.entries))
        self.assertTrue(numpy.allclose(h2.content, self.h2.content))

    def test_arithmetics(self):
        """
        """
        ha = self.h1 + self.h1
        self.assertTrue(numpy.allclose(ha.entries, 2 * self.h1.entries))
        self.assertTrue(numpy.allclose(ha.content, 2 * self.h1.content))
        self.assertTrue(numpy.allclose(ha.errors(), 2**0.5 * self.h1.errors()))
        hb = 2 * self.h1
        self.assertTrue(numpy.allclose(hb.entries, self.h1.entries))
        self.assertTrue(numpy.allclose(hb.content, ha.content))
        self.assertTrue(numpy.allclose(hb.errors(), 2 * self.h1.errors()))

    def test_scatter_plot(self):
        """
        """
        x = numpy.linspace(1, 10, 10)
        y = x ** 2
        dy = 10.
        scatter = xScatterPlot(x, y, dy, xlabel='x [a.u]', ylabel='y [a.u]')
        plt.figure('Scatter plot')
        scatter.plot()

    def test_modulation_cube(self):
        """
        """
        xedges = numpy.linspace(2., 8., 10)
        yedges = numpy.linspace(-numpy.pi, numpy.pi, 100)
        h = xModulationCube2d(xedges, yedges)
        plt.figure('Modulation cube')
        h.plot()

    def test_gpd_map(self):
        """
        """
        h = xGpdMap2d(100)
        plt.figure('GPD map')
        h.plot()

    def test_persistency(self):
        """
        """
        file_path = os.path.join(IXPEOBSSIM_DATA, 'test_histogram2d.fits')
        self.h2.save(file_path)
        hist = xHistogram2d.from_file(file_path)
        self.assertTrue(numpy.allclose(self.h2.entries, hist.entries))
        self.assertTrue(numpy.allclose(self.h2.content, hist.content))
        plt.figure('Hist from file')
        hist.plot()

    def test_kde(self):
        """
        """
        h1 = self.h1.gaussian_kde_smooth()
        plt.figure('Smoothed 1d histogram')
        h1.plot()

    def test_stats(self):
        """
        """
        sqrtn = numpy.sqrt(self.sample_size)
        delta = (self.h1.mean(0) - self.xmean) / self.xsigma * sqrtn
        self.assertTrue(abs(delta) < 5.)
        delta = (self.h2.mean(0) - self.xmean) / self.xsigma * sqrtn
        self.assertTrue(abs(delta) < 5.)
        mu = 0.5 * (self.ymin + self.ymax)
        sigma = (self.ymax - self.ymin) / numpy.sqrt(12.)
        delta = (self.h2.mean(1) - mu) / sigma * sqrtn
        self.assertTrue(abs(delta) < 5.)
        delta = self.h1.rms(0) - self.xsigma
        self.assertTrue(abs(delta) < 1.e-3)
        delta = self.h2.rms(0) - self.xsigma
        self.assertTrue(abs(delta) < 1.e-3)
        delta = self.h2.rms(1) - sigma
        self.assertTrue(abs(delta) < 1.e-3)

    def test_errors(self, size=100000):
        """
        """
        h = xHistogram1d(numpy.linspace(0, 1, 100))
        x = numpy.random.random(size)
        h.fill(x)
        self.assertTrue(numpy.allclose(h.errors(), numpy.sqrt(h.entries)))
        h = h.empty_copy()
        w = 10.
        h.fill(x, weights=w)
        self.assertTrue(numpy.allclose(h.errors(), w * numpy.sqrt(h.entries)))

    def test_find_bin(self):
        """
        """
        binning = numpy.linspace(0., 1., 3)
        h = xHistogram2d(binning, binning)
        h.content = numpy.array([[1, 2], [3, 4]])
        plt.figure('Test find bin')
        h.plot()
        self.assertEqual(h.find_bin(0.25, 0.25), (0, 0))
        self.assertEqual(h.find_bin_value(0.25, 0.25), 1)
        self.assertEqual(h.find_bin(0.25, 0.75), (0, 1))
        self.assertEqual(h.find_bin_value(0.25, 0.75), 2)
        self.assertEqual(h.find_bin(0.75, 0.25), (1, 0))
        self.assertEqual(h.find_bin_value(0.75, 0.25), 3)
        self.assertEqual(h.find_bin(0.75, 0.75), (1, 1))
        self.assertEqual(h.find_bin_value(0.75, 0.75), 4)



if __name__ == "__main__":
    unittest.main(exit=False)

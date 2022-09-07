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

'''
Module to test tdelays.py
'''

from __future__ import print_function, division

import unittest
import sys
import os

import numpy
import matplotlib.pyplot as plt

from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.srcmodel.ephemeris import xOrbitalEphemeris
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.srcmodel.tdelays import roemers_f, roemerb_f, xTDelays

if sys.flags.interactive:
    plt.ion()

# pylint: disable=invalid-name, too-many-arguments, too-many-locals

class TestTimeDelay(unittest.TestCase):

    def test_roemers_f(self):
        """
        """
        t = 0
        delay = roemers_f(t, ra=272.685, dec=-26.150)
        print("\n Delay = {}".format(delay))

    def test_roemerb_f(self):
        """
        """
        file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'parfiles', 'SAXJ1808.4-3658.par')
        ephem = xOrbitalEphemeris.from_file(file_path)
        t = 0
        delay = roemerb_f(t, ephem, ra=272.685, dec=-26.150)
        print("\n Delay = {}".format(delay))

    def test_xTDelays(self):
        """
        """
        file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'parfiles', 'SAXJ1808.4-3658.par')
        ephem = xOrbitalEphemeris.from_file(file_path)
        # Set of times
        t_values = numpy.linspace(0., 2e6, 1000)
        t = xTDelays(t_values, unit='met', name='TIME')
        t.apply_decorr(ephem, ra=45., dec=45.)
        t.apply_delay(ephem, ra=45., dec=45.)
        plt.figure('Roemer Solar System')
        Y_ = t.metvalue - t_values
        plt.plot(t_values, Y_)
        plt.title('Only earth rotation')
        plt.xlabel('MET [s]')
        plt.ylabel('Roemer Solar System Residuals [s]')

        T = xTDelays(t_values, unit='met', name='TIME')
        T.apply_decorr(ephem, ra=45., dec=45., binary=True)
        T.apply_delay(ephem, ra=45., dec=45., binary=True)
        plt.figure('Roemer Total')
        Y_ = T.metvalue - t_values
        plt.plot(t_values, Y_)
        plt.xlim((t_values.min(), t_values.max()))
        plt.title('SAXJ1808.4-3658')
        plt.xlabel('MET [s]')
        plt.ylabel('Roemer Total Residuals [s]')

        plt.figure('Roemer Binary System')
        Y_ = T.metvalue - t.metvalue
        plt.plot(t_values, Y_)
        plt.xlim((t_values.min(), t_values.max()))
        plt.title('SAXJ1808.4-3658')
        plt.xlabel('MET [s]')
        plt.ylabel('Roemer Binary System Residuals [s]')

if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

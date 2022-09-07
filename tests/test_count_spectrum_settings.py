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

from __future__ import print_function, division


"""Unit test for the count spectrum facility.
"""

import numpy
import unittest
import sys

from ixpeobssim.irf import load_arf
from ixpeobssim.srcmodel.spectrum import xCountSpectrum
from ixpeobssim.utils.matplotlib_ import plt
import ixpeobssim.config.toy_livetime_time as input_model

if sys.flags.interactive:
    plt.ion()



class TestCountSpectrum(unittest.TestCase):

    """Unit test for xCountSpectrum.
    """

    def test_toy_livetime_time(self, slice_energy=3., num_points=500):
        """Create a count spectrum
        """
        src = input_model.src
        print(src)
        aeff = load_arf(du_id=1)
        start_met = input_model.START_MET
        stop_met = input_model.START_MET + input_model.DURATION
        kwargs = dict(start_met=start_met, duration=input_model.DURATION)
        # Default settings...
        time_grid = src.sampling_time_grid(*src.parse_time_kwargs(**kwargs))
        fine_time_grid = numpy.linspace(start_met, stop_met, 1000)
        count_spectrum = src.create_count_spectrum(aeff, time_grid, **kwargs)
        plt.figure('Default count spectrum')
        count_spectrum.plot()
        plt.figure('count spectrum slice')
        plt.plot(fine_time_grid, count_spectrum(slice_energy, fine_time_grid))
        # Modified settings.
        src.set_count_spectrum_params(num_points, 3, 1)
        time_grid = src.sampling_time_grid(*src.parse_time_kwargs(**kwargs))
        count_spectrum = src.create_count_spectrum(aeff, time_grid, **kwargs)
        plt.figure('Modified count spectrum')
        count_spectrum.plot()
        plt.figure('count spectrum slice')
        self.assertEqual(count_spectrum.y.shape, (num_points,))
        plt.plot(fine_time_grid, count_spectrum(slice_energy, fine_time_grid))



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

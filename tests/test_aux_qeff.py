#!/usr/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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

from ixpeobssim.irfgen.auxiliary import AUX_VERSION, AUX_REFERENCE_PRESSURE,\
    AUX_ABSORPTION_LABELS, AUX_WEIGHT_NAMES, lines_qeff_file_path,\
    allx_qeff_file_path
from ixpeobssim.irfgen.auxiliary import load_qeff_table
from ixpeobssim.irfgen.gpd import xQeffDataInterface
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color

if sys.flags.interactive:
    plt.ion()


class TestAuxQeff(unittest.TestCase):

    """Test the quantum efficieny auxiliary files.
    """

    @classmethod
    def setUpClass(cls, aux_version=AUX_VERSION):
        """Load the effective area data for all the relevant weights
        """
        cls.qeff_data = {}
        cls.line_data = {}
        cls.allx_data = {}
        for weight_name in AUX_WEIGHT_NAMES:
            cls.qeff_data[weight_name] = xQeffDataInterface(weight_name, aux_version)
            cls.line_data[weight_name] = load_qeff_table(lines_qeff_file_path(weight_name, aux_version))
            cls.allx_data[weight_name] = load_qeff_table(allx_qeff_file_path(weight_name, aux_version))

    def _test_raw_data_base(self, weight_name):
        """Plot the raw data (in the line and allx flavors).
        """
        plt.figure('Raw qeff data (weights = %s)' % weight_name)
        line_data = self.line_data[weight_name]
        allx_data = self.allx_data[weight_name]
        for label in AUX_ABSORPTION_LABELS:
            y = line_data['QEFF_%s' % label.upper()]
            dy = line_data['QEFF_%s_ERR' % label.upper()]
            plt.errorbar(line_data['ENERGY'], y, dy, fmt='o')
            y = allx_data['QEFF_%s' % label.upper()]
            plt.plot(allx_data['ENERGY'], y, color=last_line_color(), label=label)
        setup_gca(xlabel='Energy [keV]', ylabel='Quantum efficiency @ %d mbar' % AUX_REFERENCE_PRESSURE,
                  xmin=1., ymin=1.e-5, logy=True, grids=True, legend=True)

    def test_raw_data(self):
        """Plot all the raw data.
        """
        for weight_name in AUX_WEIGHT_NAMES:
            self._test_raw_data_base(weight_name)

    def _test_interface_residuals_base(self, weight_name):
        """
        """
        qeff_data = self.qeff_data[weight_name]
        line_data = self.line_data[weight_name]
        allx_data = self.allx_data[weight_name]

        fmt = dict(xlabel='Energy [keV]', grids=True, xmin=1., ymin=-0.1, ymax=0.1,
                   ylabel='Fractional residuals @ %d mbar' % AUX_REFERENCE_PRESSURE)
        plt.figure('qeff window residuals (weights = %s)' % weight_name)
        x = allx_data['ENERGY']
        y = allx_data['QEFF_WIN']
        dy = allx_data['QEFF_WIN_ERR']
        model = qeff_data.window_quantum_efficiency(x)
        res = y - model
        res /= model
        dy /= model
        plt.errorbar(x, res, dy, fmt='o')
        setup_gca(**fmt)

        plt.figure('qeff DME residuals (weights = %s)' % weight_name)
        x = allx_data['ENERGY']
        y = allx_data['QEFF_DME']
        dy = allx_data['QEFF_DME_ERR']
        model = qeff_data.dme_quantum_efficiency(x, pressure=AUX_REFERENCE_PRESSURE, contaminants=None)
        res = y - model
        res /= model
        dy /= model
        plt.errorbar(x, res, dy, fmt='o')
        setup_gca(**fmt)

        plt.figure('qeff GEM residuals (weights = %s)' % weight_name)
        x = allx_data['ENERGY']
        y = allx_data['QEFF_GEM']
        dy = allx_data['QEFF_GEM_ERR']
        model = qeff_data.gem_quantum_efficiency(x, pressure=AUX_REFERENCE_PRESSURE, contaminants=None)
        res = y - model
        res /= model
        dy /= model
        plt.errorbar(x, res, dy, fmt='o')
        setup_gca(**fmt)

    def test_interface_residuals(self):
        """
        """
        for weight_name in AUX_WEIGHT_NAMES:
            self._test_interface_residuals_base(weight_name)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

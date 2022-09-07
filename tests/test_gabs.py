#!/urs/bin/env python
#
# Copyright (C) 2015--2016, the ixpeobssim team.
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


import sys
import os
import unittest
import numpy

from ixpeobssim.srcmodel.gabs import mapped_column_density_HI,\
    xInterstellarAbsorptionModel
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim import IXPEOBSSIM_SRCMODEL, IXPEOBSSIM_DOC_FIGURES

if sys.flags.interactive:
    plt.ion()


class TestInterstellarAbsorption(unittest.TestCase):

    """Unit test for the interstellar absorption.
    """

    def test_transmission(self):
        """
        """
        model = xInterstellarAbsorptionModel()
        for nh, label in [(1e21, '1e21'), (1e22, '1e22'), (1e23, '1e23')]:
            trans = model.transmission_factor(nh)
            file_name = 'xspec_wabs_%s.txt' % label
            file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'ascii', file_name)
            elo, ehi, txspec = numpy.loadtxt(file_path, unpack=True)
            energy = numpy.sqrt(elo*ehi)
            _mask = (txspec > 1.e-2)
            energy = energy[_mask]
            txspec = txspec[_mask]
            tixpeobssim = trans(energy)
            delta = abs(txspec - tixpeobssim)/txspec
            delta_ave = numpy.average(delta)
            self.assertTrue(delta_ave < 0.01, 'average difference %s' % delta_ave)

    def test_plots(self):
        """
        """
        model = xInterstellarAbsorptionModel()
        plt.figure()
        model.xsection.plot(logx=True, logy=True)

        plt.figure()
        model.xsection_ecube().plot(logx=True)

        plt.figure()
        ra, dec = 10.684, 41.269
        column_density = mapped_column_density_HI(ra, dec, 'LAB')
        trans = model.transmission_factor(column_density)
        trans.plot(logx=True, label='$n_H = $%.3e' % column_density)
        plt.legend(loc='upper left')

        plt.figure()
        for column_density in [1.e20, 1.e21, 1.e22, 1.e23]:
            trans = model.transmission_factor(column_density)
            trans.plot(logx=True, label='$n_H = $%.1e' % column_density)
        plt.legend(loc='upper left')



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

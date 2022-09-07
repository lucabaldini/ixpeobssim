#!/usr/bin/env python
#
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


import unittest
import os

import numpy

import ixpeobssim.core.pipeline as pipeline
from ixpeobssim import IXPEOBSSIM_CONFIG_REG
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.astro import read_ds9
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt


class TestSelect(unittest.TestCase):

    """Unit test for simulation and analysis pipeline(s).
    """

    @classmethod
    def setUpClass(cls):
        """Run a short Cas A simulation.

        We want an extended source to test ds9 region files, but we also fold it
        in phase with arbitrary ephemeris to test the phase selection.
        """
        pipeline.reset('toy_casa', overwrite=True)
        file_list = pipeline.xpobssim(duration=1000., saa=False, occult=False)
        cls._file_list = pipeline.xpphase(*file_list, nu0=1.)
        evt_file = xEventFile(file_list[0])
        cls._start_met = evt_file.start_met()
        cls._stop_met = evt_file.stop_met()

    def _test_base(self, kwargs, invert_switch):
        """This is checking that the the sum of the events in the complementary
        selections amount to the original total.
        """
        kwargs1 = kwargs
        kwargs2 = kwargs1.copy()
        kwargs2[invert_switch] = True
        kwargs2['suffix'] = '%sinv' % kwargs1['suffix']
        file_list1 = pipeline.xpselect(*self._file_list, **kwargs1)
        file_list2 = pipeline.xpselect(*self._file_list, **kwargs2)
        for file_path, file_path1, file_path2 in zip(self._file_list, file_list1, file_list2):
            f = xEventFile(file_path)
            f1 = xEventFile(file_path1)
            f2 = xEventFile(file_path2)
            n = f.num_events()
            n1 = f1.num_events()
            n2 = f2.num_events()
            logger.info('Event splitting: %d -> %d + %d', n, n1, n2)
            self.assertEqual(n, n1 + n2)

    def test_time(self):
        """Test the time selection.
        """
        kwargs = dict(tmax=0.5 * (self._start_met + self._stop_met), suffix='tsel')
        self._test_base(kwargs, 'tinvert')

    def test_phase(self):
        """Test the phase selection.
        """
        kwargs = dict(phasemin=0.5, suffix='phasesel')
        self._test_base(kwargs, 'phaseinvert')

    def test_energy(self):
        """Test the energy selection.
        """
        kwargs = dict(emin=2, emax=6, suffix='esel')
        self._test_base(kwargs, 'einvert')

    def test_reg(self):
        """Test the region selection.
        """
        reg_file_path = os.path.join(IXPEOBSSIM_CONFIG_REG, 'casa_selected.reg')
        kwargs = dict(regfile=reg_file_path, suffix='regsel')
        self._test_base(kwargs, 'reginvert')




if __name__ == '__main__':
    unittest.main()

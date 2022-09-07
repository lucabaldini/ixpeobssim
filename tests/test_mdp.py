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


"""Unit test for the MDP facilities.
"""

import unittest

import numpy

from ixpeobssim.evt.mdp import mdp99, xMDPRecord, xMDPTable
from ixpeobssim.utils.misc import pairwise


class TestMDP(unittest.TestCase):

    """Unit test for the MDP facilities.
    """

    def test_mdp99(self):
        """
        """
        self.assertAlmostEqual(mdp99(0.5, 750000), 0.009911949421447495)

    def test_print(self):
        """
        """
        emin=2.
        emax=8.
        print(xMDPRecord.empty(emin, emax))
        print(xMDPRecord(emin, emax, 0.3, 100000))
        print(xMDPRecord(emin, emax, 0.3, 100000, 2000))

    def test_record_addition(self):
        """
        """
        emin=2.
        emax=8.
        record0 = xMDPRecord.empty(emin, emax)
        record1 = xMDPRecord(emin, emax, 0.3, 100000, 0)
        record2 = record1 + record1
        print(record0)
        print(record1)
        print(record2)
        self.assertAlmostEqual(record1.mdp, record2.mdp * numpy.sqrt(2.))
        self.assertAlmostEqual((record1 + record0).mdp, record1.mdp)

    def test_table_addition(self):
        """
        """
        ebinning = numpy.linspace(2., 8., 5)
        duration = 10000.
        table0 = xMDPTable.empty(duration, ebinning)
        table1 = xMDPTable(duration)
        for emin, emax in pairwise(ebinning):
            table1.add_row(emin, emax, 0.3, 100000)
        print(table0)
        print(table1)
        print(table0 + table1)
        self.assertTrue(table1.total_counts() == 400000)

    def test_scale(self):
        """
        """
        emin=2.
        emax=8.
        record = xMDPRecord(emin, emax, 0.3, 100000, 0)
        print(record)
        mdp = record.mdp
        record.scale(2.)
        print(record)
        self.assertAlmostEqual(record.mdp, mdp / numpy.sqrt(2.))



if __name__ == '__main__':
    unittest.main()

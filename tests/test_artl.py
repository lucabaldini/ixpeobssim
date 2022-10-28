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

"""Unit tests for the targets.artl module.
"""


from __future__ import print_function, division

import sys
import unittest

from ixpeobssim.targets.__artl__ import SOURCE_LEGEND_HANDLES
from ixpeobssim.targets.artl import TARGET_DICT, xObservation
from ixpeobssim.utils.matplotlib_ import plt

if sys.flags.interactive:
    plt.ion()


class TestArtl(unittest.TestCase):

    """Unit test for targets.artl module.
    """

    def test_targets(self):
        """
        """
        fig = plt.figure('Target positions', figsize=(12., 8.))
        ax = fig.add_subplot(111, projection='aitoff')
        for target in TARGET_DICT.values():
            print(target)
            target.plot_galactic_pos()
        ax.grid(True)
        plt.gca().legend(handles=SOURCE_LEGEND_HANDLES, loc=(-0.15, 0.91))

    def test_observation(self):
        """
        """
        obs = xObservation('Cas A', 11.57, 1001301, '2022-01-11T11:45', '2022-01-29T12:34', '')
        print(obs)
        obs = xObservation('Mrk 501', 1.16, 2004501, '2022-12-25T', '2022-12-27T', '')
        print(obs)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

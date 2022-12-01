#!/usr/bin/env python
#
# Copyright (C) 2018, the ixpeobssim team.
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
import numpy

from ixpeobssim.evt.clustering import region_query_factory, PixelStatus, DBscan
from ixpeobssim.evt.display import xL1EventFile, xL1Event
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt


@unittest.skip('This has a hard-coded path')
class TestClocking(unittest.TestCase):

    """Unit test for the clustering.
    """
    @classmethod
    def setUpClass(cls, du_id=1):
        """Setup---load the event file.
        """
        file_path = '/media/alberto/TOSHIBA EXT/xpe/xpedata/3C_279/01005701/event_l1/ixpe01005701_det1_evt1_v03.fits'
        cls.evt_file = xL1EventFile(file_path)

    def test_region_query(self, threshold=18):
        """
        """
        for i in range(1, 10):
            evt = self.evt_file[i]
            # Pick a pixel from the middle of the event
            idx = int(evt.size/2)
            # Get its offset coordinates
            cols, rows = evt.serial_readout_coordinates()
            col, row = cols[idx], rows[idx]
            # Switch from absolute (i.e. referred to the entre CHIP) to relative
            # (i.e. referred to the current ROI) offset coordinates
            rel_col = col - evt.min_col
            rel_row = row - evt.min_row
            # Create the region query function for this event
            region_query = region_query_factory(evt)
            # Call it on the pxiel
            neighbors = region_query(idx, threshold)
            # Check that the number of neighbors is smaller or equal than 6
            self.assertTrue(len(neighbors) <= 6)
            # Check that all the neighbors found are over threshold
            for neighbor in neighbors:
                _col, _row = cols[neighbor], rows[neighbor]
                _col -= evt.min_col
                _row -= evt.min_row
                self.assertTrue(evt.pha[_row, _col] >= threshold)

    def test_clustering(self, threshold=20, min_density_points=4,
                        min_cluster_size=6):
        """
        """
        # Initialize the clustering facility
        clustering = DBscan(threshold, min_density_points, min_cluster_size)
        for i in range(1, 10):
            evt = self.evt_file[i]
            # Initialize the region query function for the current event
            region_query = region_query_factory(evt)
            # Get the pha array from the event - we need it flat for the clustering
            input_data = evt.pha.flatten()
            # Initialize the output array - this can actually be initialized
            # with any value, as it is going to be filled by the clustering loop
            output_ids = numpy.full(input_data.shape, -999)
            # Run the clustering on the current event
            clusters_info = clustering.run(input_data, output_ids, region_query)
            # Check that all the pixels have been visited by the algorithm
            self.assertTrue(numpy.all(output_ids > PixelStatus.UNDEFINED))
            clusters_total_pha = []
            for i, cluster in enumerate(clusters_info):
                # Check that cluster numbering is in increasing order, no gaps
                self.assertEqual(i, cluster.cluster_id)
                # Check that the cluster size is above the minimum size
                self.assertTrue(cluster.num_pixels >= min_cluster_size)
                # Check that all the pixels in the cluster are above threshold
                _mask = output_ids == cluster.cluster_id
                self.assertTrue(numpy.all(input_data[_mask] >= threshold))
                # Check that the total PHA of the pixels match the number
                # written in the ClusterInfo object
                total_pha = numpy.sum(input_data[_mask])
                self.assertEqual(total_pha, cluster.pulse_height)
                clusters_total_pha.append(total_pha)
            # Check that the clusters are indeed sorted by PHA
            self.assertEqual(sorted(clusters_total_pha, reverse=True),
                             clusters_total_pha)



if __name__ == '__main__':
    unittest.main()

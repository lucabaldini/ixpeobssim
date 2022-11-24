#/***********************************************************************
#Copyright (C) 2022 the GPD team.
#
#For the license terms see the file LICENSE, distributed along with this
#software.
#
#This program is free software; you can redistribute it and/or modify it
#under the terms of the GNU General Public License as published by the
#Free Software Foundation; either version 2 of the License, or (at your
#option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#General Public License for more details.
#
#You should have received a copy of the GNU General Public License along
#with this program; if not, write to the Free Software Foundation Inc.,
#51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#***********************************************************************/
from __future__ import annotations

import numpy
import enum
from dataclasses import dataclass


"""
"""


@dataclass
class CubeCoordinate:
    """ Cube coordinate
    """

    q: int
    r: int
    s: int

    def __add__(self, other):
        """+= operator."""
        return CubeCoordinate(self.q + other.q,
                              self.r + other.r,
                              self.s + other.s)

    def __sub__(self, other):
        """+= operator."""
        return CubeCoordinate(self.q - other.q,
                              self.r - other.r,
                              self.s - other.s)


class EvenRTransforms:
    """Class implementing transformation from offset to cubic coordinates
    an vice versa. This depends essentially from the ASIC grid layout.
    Coordinates transformation are largely mutuated from
    http://www.redblobgames.com/grids/hexagons/implementation.htm
    """

    @staticmethod
    def cubic_to_offset(c):
        """ Cube to Offset transformation for Even R grids.

        Arguments
        ---------
        c : CubeCoordinate
            The input cubic coordinate.
        """
        col = c.q + int((c.r + (c.r & 1)) / 2)
        row = c.r
        return col, row

    @staticmethod
    def offset_to_cubic(col, row):
        """ Offset to Cube transformation for Even R grids.

        Arguments
        ---------
        col: int
            Input column (in CHIP coordinates system).

        row: int
            Input row (in CHIP coordinates system).
        """
        q = col - int((row + (row & 1)) / 2)
        r = row
        s = -q - r
        return CubeCoordinate(q, r, s)


def find_neighbors_cubic(coords):
    """ Get a list of the cubic coordinates of the neighbors of
    a given point.

    Arguments
    ---------
    coords: CubeCoordinate
        The cubic coordinate of the pixel.
    """
    return [coords + CubeCoordinate(+1, -1,  0),
            coords + CubeCoordinate(+1,  0, -1),
            coords + CubeCoordinate( 0, +1, -1),
            coords + CubeCoordinate(-1, +1,  0),
            coords + CubeCoordinate(-1,  0, +1),
            coords + CubeCoordinate( 0, -1, +1)]


def find_neighbors(col, row):
    """ Get a list of the offset coordinates of the neighbors of a given
    point.

    Arguments
    ---------
    col: int
        Input pixel column (in CHIP coordinates system).

    row: int
        Input pixel row (in CHIP coordinates system).
    """
    cc = EvenRTransforms.offset_to_cubic(col, row)
    return [EvenRTransforms.cubic_to_offset(id) \
                for id in find_neighbors_cubic(cc)]


def region_query_factory(event):
    """ Factory for the auxiliary region_query function used by the
    clustering. We use the closure to cache some data from the current event,
    namely the two arrays of offset coordinates (in serial readout order)
    and the corresponding 2D array of serial readout indices.

    Arguments
    ---------
    event: xL1Event
        The current event.
    """
    cols, rows = event.serial_readout_coordinates()
    indices = event.serial_readout_indices()

    def region_query(index, threshold):
        """ Find the over threshold neighbors of a given pixel in the current
        event.

        Arguments
        ---------
        index: int
            The pixel readout index.

        threshold: float
            The clustering threshold.
        """
        neighbor_offset_coords = find_neighbors(cols[index], rows[index])
        neighbors = []
        for col, row in neighbor_offset_coords:
            if not event.coordinates_in_roi(col, row):
                continue
            # Switch from absolute (i.e. referred to the entre CHIP) to relative
            # (i.e. referred to the current ROI) offset coordinates
            col -= event.min_col
            row -= event.min_row
            if event.pha[row, col] > threshold:
                neighbors.append(indices[row, col])
        return neighbors

    return region_query


class PixelStatus(enum.IntEnum):
    """Enum for the different status of a point (pixel) while running the
    DBSCAN algorithm. See the DBscan class for a detailed explanation.
    Note that we need all the corresponding values to be negative, as we
    reserve positive values for the cluster IDs."""
    UNDEFINED = -3 # Pixel yet to be visited by the clustering algorithm
    SUBTHRESHOLD = -2 # Pixel not assigned to any cluster (under threshold)
    NOISE = -1 # Pixel not assigned to any cluster (over threshold but isolated)


@dataclass
class ClusterInfo:
    """ Class for storing information about a cluster."""
    cluster_id : int
    num_pixels : int = 0
    pulse_height : float = 0.

    def update(self, pha):
        """ Update the cluster information with a new pixel.

        Arguments
        ---------
        pixel_value:  float
            The pixel PHA.
        """
        self.num_pixels += 1
        self.pulse_height += pha

    def __gt__(self, other):
        """ Define the rules for cluster sorting (in order of priority):
        1) greater pulse invariant
        2) greater number of pixels
        """
        if self.pulse_height < other.pulse_height:
            return False
        if self.num_pixels < other.num_pixels:
            return False
        return True

    def __str__(self):
        return f'Cluster id {self.cluster_id}: {self.num_pixels} pixels, '\
               f'{self.pulse_height} ADC counts total'


class DBscan:
    """Class implmenting the DBSCAN clustering algorithm
    See https://en.wikipedia.org/wiki/DBSCAN for an introduction to DBSCAN

    Original work:
    Ester, Martin; Kriegel, Hans-Peter; Sander, Jörg; Xu, Xiaowei (1996)
    "A density-based algorithm for discovering clusters in large spatial
    databases with noise." Proceedings of the Second International Conference
    on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231.
    CiteSeerX 10.1.1.121.9220 Freely accessible. ISBN 1-57735-004-9.

    The DBSCAN algorithm is based on the dfinition of a neighborood for each
    point. The points are classified as core points, (density-) reachable
    points and outliers, as follows:

    - A point p is a core point if at least min_denisty points are in its
      neighborood (including p).
    - A point q is directly reachable from a core point p if it is in the
      neighborood of p. Points are only said to be directly reachable
      from core points.
    - A point q is reachable from (or density connected to) p if there is a
      path connecting p and q passing only through points each directly
      reachable from the previous one.
      Note that this implies that the initial point and all points on the
      path must be core points, with the possible exception of q.

    All points not reachable from any other point are outliers or noise points.

    In the original algorithm, the neighborood of a point is defined by
    specyifing a radius around the point. In our case, however, a pixel can
    only be neighbor with its ajacent pixels on the hexagonal grid, and we
    replace the radius parameter with a threshold for noise suppression.
    So the neighborood of a pixel is defined by the number of hexagonally
    adjacent pixels which are above the noise suppression threshold.

    Arguments
    ---------
    threshold : float
        The clustering threshold.

    min_density_points : int
        The minimum number of over-threshold neighbors required for a pixel to
        be classifed as CORE (including the pixel itself).

    min_cluster_size: int
        The minimum size of a cluster.
    """
    def __init__(self, threshold, min_density_points=4, min_cluster_size=6):
        """ Class constructor.
        """
        self.threshold = threshold
        self.min_density_points = min_density_points
        self.min_cluster_size = min_cluster_size

    def run(self, input_data, output_ids, region_query):
        """ The clustering workhorse function. Fill the output_data array with
        the cluster id for each pixel and also returns a sorted list of
        ClusterInfo objects.

        Arguments
        ---------
        input_data : numpy.array
            The 1-dimensional pha values (they would typically be in serial
            readout order, altough that is not required, as long as the
            index is consistent with the logic of the region_query function)

        output_ids : numpy.array
            A 1-dimensional array of the same shape of the input_data, which
            will be filled with the cluster id for each pixel.

        region_query: func (int, float) -> List[int]
            A function returning the list of indices of the over-threshold
            neighbors of a pixel in the current event given the pixel index
            and a threshold

        """
        cluster_ids = numpy.full(output_ids.shape, PixelStatus.UNDEFINED)
        # Create an empty list of ClusterInfo objects
        clusters  = []
        # Main loop
        for pixel_id, pixel_value in enumerate(input_data):
            # Skip already visited pixels
            if cluster_ids[pixel_id] != PixelStatus.UNDEFINED:
                continue
            # Flag under threshold pixels
            # TODO: we may consider accepting under threshold pixels, as long as
            # they are surrounded by pixels in the cluster
            if pixel_value <= self.threshold:
                cluster_ids[pixel_id] = PixelStatus.SUBTHRESHOLD
                continue
            # Get an array of adjacent, over threshold pixels
            neighbors = region_query(pixel_id, self.threshold)
            # Check wheter the pixel is isolated
            if len(neighbors) < self.min_density_points:
                # Flag the pixel as NOISE (over threshold, but isolated)
                cluster_ids[pixel_id] = PixelStatus.NOISE
                continue
            # If not, start building a new cluster from this pixel
            # Create the new cluster (with incremental id)
            cluster = ClusterInfo(cluster_id=len(clusters))
            # Add the pixel to the cluster
            cluster_ids[pixel_id] = cluster.cluster_id
            cluster.update(pixel_value)
            # Loop over adiacent pixels (we can't use an iterator because of the
            # array being expanded inside the loop, which would invalidate it)
            i = 0
            while i < len(neighbors):
                # Get the id of the currently visted pixel
                _id = neighbors[i]
                i += 1
                # If the pixel has already been flagged as noise, now merge it
                # into to the cluster - but do not expand from it (it is a
                # border pixel)
                if cluster_ids[_id] == PixelStatus.NOISE:
                    cluster_ids[_id] = cluster.cluster_id
                    cluster.update(input_data[_id])
                    continue
                # If the pixel status is not undefined skip it
                if cluster_ids[_id] != PixelStatus.UNDEFINED:
                    continue
                # Assign the pixel to the current cluster
                cluster_ids[_id] = cluster.cluster_id
                cluster.update(input_data[_id])
                # Get the array of its neighbors
                _neighbors = region_query(_id, self.threshold)
                # If the density condition is satisfied, expand the list of
                # candidate pixels to visit with the neighbors of the current
                # pixel
                if len(_neighbors) >= self.min_density_points:
                    for pixel in _neighbors:
                        if pixel not in neighbors:
                            neighbors.append(pixel)
            # Finally add the new cluster to the array of clusters found
            clusters.append(cluster)
        # Now get rid of all the clusters below the minimum size allowed,
        # marking their pixels as NOISE
        good_clusters = []
        for cluster in clusters:
            if cluster.num_pixels < self.min_cluster_size:
                cluster_mask = (cluster_ids == cluster.cluster_id)
                cluster_ids[cluster_mask] = PixelStatus.NOISE
            else:
                good_clusters.append(cluster)
        # Now sort clusters in descending order (so that the first cluster has
        # the highest pulse invariant)
        good_clusters.sort(reverse=True)
        # Finally set the indices in the output array to the corresponding
        # sorted cluster id
        numpy.copyto(output_ids, cluster_ids)
        for idx, cluster in enumerate(good_clusters):
            cluster_mask = (cluster_ids == cluster.cluster_id)
            output_ids[cluster_mask] = idx
            # Update the cluster id with the sorted number
            cluster.cluster_id = idx
        return good_clusters

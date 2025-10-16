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

"""Data structures related to the good time intervals.
"""

from __future__ import print_function, division

from enum import IntEnum

import numpy

from ixpeobssim.utils.logging_ import logger

# pylint: disable=invalid-name


class EdgeType(IntEnum):
    """
    Enumeration representing the types of edges that can occur in Good
    Time Intervals (GTIs) and Bad Time Intervals (BTIs).

    Members:
        GTI_START: Beginning of a good time interval.
        GTI_STOP: End of a good time interval.
        BTI_START: Beginning of a bad time interval.
        BTI_STOP: End of a bad time interval.
    """
    GTI_START = 0
    GTI_STOP = 1
    BTI_START = 2
    BTI_STOP = 3


class UnexpectedEdgeType(RuntimeError):
    """
    Custom exception raised when an unexpected sequence of GTI/BTI edges is
    encountered (e.g. two GTI start ina row without a GTI stop in between).

    Attributes:
        msg (str): A detailed error message describing the unexpected edge and
        current state.
    """
    def __init__(self, edge_time, edge_type, in_old_gti, in_bti):
        """
        Initialize the exception with the edge time, type, and current interval
        state.

        Args:
            edge_time: The time of the edge event.
            edge_type: The type of the edge (GTI_START, GTI_STOP, etc.).
            in_old_gti (bool): True if currently inside an old GTI.
            in_bti (bool): True if currently inside a BTI.
        """
        self.msg = f'Unexpected edge type {str(EdgeType(edge_type))} at time '\
                   f'{edge_time}, when "in_old_gti" is {in_old_gti} and '\
                   f'"in_bti" is {in_bti}'

    def __str__(self):
        return self.msg


class xGTIList(list):

    """Small convenience class representing a list of good time interval.

    The interfaces are fairly minimal,  but since we use this in quite
    different places, it was handy to collect the useful stuff in one place.
    """

    def __init__(self, start_met, stop_met, *gtis):
        """
        """
        list.__init__(self)
        self.start_met = start_met
        self.stop_met = stop_met
        self.span = self.stop_met - self.start_met
        for start, stop in gtis:
            self.append_gti(start, stop)

    def append_gti(self, start, stop):
        """Append a new GTI to the list.
        """
        assert start >= self.start_met
        assert stop <= self.stop_met
        self.append((start, stop))

    def total_good_time(self):
        """Return the total good time.
        """
        return sum(stop - start for start, stop in self)

    def start_mets(self):
        """Return the value of the start MET values for the GTI.
        """
        return [start for start, stop in self]

    def stop_mets(self):
        """Return the value of the stop MET values for the GTI.
        """
        return [stop for start, stop in self]

    def all_mets(self):
        """Return all the MET values corresponding to the start or stop of the
        GTIs in the list.

        Note that Python supports the call to sum() with start set by keyword
        argument only since version 3.8, so we are refraining from that, here.
        """
        return sum([[start, stop] for start, stop in self], [])

    def complement(self):
        """Return the logical complement of the GTI list in the relevant time
        span, i.e., the list of time intervals between the start and the end
        of observation that are *not* good time intervals.

        These are the ones where we do not take celestial data and we can, e.g.,
        take calibration data.

        See https://stackoverflow.com/questions/16789776/ for the use of the
        iterator over the list.
        """
        # This is an iterator to all the bounds of the GTIs except for the
        # start of the first GTI and the end of the last.
        iterator = iter(self.all_mets()[1:-1])
        return [(start, next(iterator)) for start in iterator]

    def filter_event_times(self, time_):
        """Filter a given array of event times and return a reduced array
        only containing the times within the good time intervals in the list,
        along with the corresponding boolean mask. The latter, in turn, can be
        used to filter ancillary related columns, e.g., the phase in simulations
        of periodic sources.
        """
        logger.info('Filtering %d event times according to the GTIs...', len(time_))
        mask = numpy.zeros(time_.shape, dtype=bool)
        for (start, stop) in self:
            mask[numpy.logical_and(time_ >= start, time_ <= stop)] = True
        time_ = time_[mask]
        logger.info('Done, %d entries remaining.', len(time_))
        return time_, mask

    def update(self, bti_start, bti_stop):
        """
        Update Good Time Intervals (GTIs) by removing Bad Time Intervals (BTIs)
        from them.
 
        This function merges and processes the GTI and BTI edges in chronological
        order, producing new GTIs that exclude the periods defined as BTIs.
 
        Args:
            bti_start (numpy.ndarray): Array of start times for BTIs.
            bti_stop (numpy.ndarray): Array of stop times for BTIs.
 
        Returns:
            Tuple[numpy.ndarray, numpy.ndarray]: Arrays of start and stop times for
            the new GTIs.
        """
        gti_start = numpy.array(self.start_mets())
        gti_stop = numpy.array(self.stop_mets())
 
        # Combine all edge times and types
        edge_times = numpy.hstack((gti_start, gti_stop, bti_start, bti_stop))
        edge_types = numpy.array(
            [EdgeType.GTI_START] * len(gti_start) +
            [EdgeType.GTI_STOP] * len(gti_stop) +
            [EdgeType.BTI_START] * len(bti_start) +
            [EdgeType.BTI_STOP] * len(bti_stop)
        )

        # Sort edges by time (stable sort ensures start precedes stop if times are equal)
        idx = numpy.argsort(edge_times, kind='mergesort')
        edge_times = edge_times[idx]
        edge_types = edge_types[idx]

        # State trackers
        in_old_gti = False
        in_bti = False
        new_gti_start = []
        new_gti_stop = []

        # Process each edge in time order
        for edge_time, edge_type in zip(edge_times, edge_types):
            # Case 1: outside both GTI and BTI
            if not in_old_gti and not in_bti:
                if edge_type == EdgeType.GTI_START:
                    new_gti_start.append(edge_time)
                    in_old_gti = True
                elif edge_type == EdgeType.BTI_START:
                    in_bti = True
                else:
                    raise UnexpectedEdgeType(edge_time, edge_type, in_old_gti,
                                             in_bti)
            # Case 2: inside BTI but not GTI
            elif not in_old_gti and in_bti:
                if edge_type == EdgeType.GTI_START:
                    in_old_gti = True
                elif edge_type == EdgeType.BTI_STOP:
                    in_bti = False
                else:
                    raise UnexpectedEdgeType(edge_time, edge_type, in_old_gti,
                                             in_bti)
            # Case 3: inside GTI but not BTI
            elif in_old_gti and not in_bti:
                if edge_type == EdgeType.GTI_STOP:
                    in_old_gti = False
                    new_gti_stop.append(edge_time)
                elif edge_type == EdgeType.BTI_START:
                    in_bti = True
                    new_gti_stop.append(edge_time)  # End current GTI at BTI start
                else:
                    raise UnexpectedEdgeType(edge_time, edge_type, in_old_gti,
                                             in_bti)
            # Case 4: inside both GTI and BTI
            else:
                if edge_type == EdgeType.GTI_STOP:
                    in_old_gti = False
                elif edge_type == EdgeType.BTI_STOP:
                    in_bti = False
                    new_gti_start.append(edge_time)  # Start new GTI after BTI ends
                else:
                    raise UnexpectedEdgeType(edge_time, edge_type, in_old_gti,
                                             in_bti)

        return numpy.array(new_gti_start), numpy.array(new_gti_stop)
 
    def __str__(self):
        """String formatting.
        """
        text = 'List of good time intervals (%.3f s total over %.3f s span):' %\
            (self.total_good_time(), self.span)
        for i, (start, stop) in enumerate(self):
            text = '%s\n[%3d] (%.3f--%.3f) or (%.3f--%.3f)' %\
                (text, i + 1, start, stop, start - self.start_met,
                 stop - self.start_met)
        return text


class xSimpleGTIList(xGTIList):

    """Subclass of xGTIList with a single GTI.
    """

    def __init__(self, start_met, stop_met):
        """Constructor.
        """
        xGTIList.__init__(self, start_met, stop_met, (start_met, stop_met))



class xUberGTIList(xGTIList):

    """Subclass of xGTIList with a single GTI including all the MET from -inf to inf.
    """

    def __init__(self):
        """Constructor.
        """
        start_met = -numpy.inf
        stop_met = numpy.inf
        xGTIList.__init__(self, start_met, stop_met, (start_met, stop_met))

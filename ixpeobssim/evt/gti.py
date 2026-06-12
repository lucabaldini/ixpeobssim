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



class UnexpectedEdgeType(RuntimeError):
    """
    Custom exception raised when an unexpected sequence of GTI edges is
    encountered (e.g. two GTI start in a row without a GTI stop in between).
    """
    pass


class xGTIList(list):

    """Class representing a list of good time interval.

    The interfaces are fairly minimal, but since we use this in quite
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

    @classmethod
    def from_arrays(cls, start_met, stop_met, tstarts, tstops):
        """
        Initialize from two arrays (of start and stop).
        This is essentially converting the two arrays of start and stop to 
        an array of tuples (start, stop), as required by the constructor of
        the class.
        """
        gtis = []
        for start, stop in zip(tstarts, tstops):
            gtis.append([start, stop])
        return cls(start_met, stop_met, *gtis)

    @staticmethod
    def validate_intervals(start, stop):
        if len(start) != len(stop):
            raise ValueError('start/stop length mismatch')
        if not numpy.all(start < stop):
            raise ValueError('non-positive interval detected')

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

    def gti_mask(self, time_):
        """
        Return a boolean mask selecting the times falling inside the GTI in the 
        list.
        
        Args
        ---------

        time_ : array of times to select
        """
        mask = numpy.zeros(time_.shape, dtype=bool)
        for (start, stop) in self:
            mask[numpy.logical_and(time_ >= start, time_ <= stop)] = True
        return mask

    def filter_event_times(self, time_):
        """Filter a given array of event times and return a reduced array
        only containing the times within the good time intervals in the list,
        along with the corresponding boolean mask. The latter, in turn, can be
        used to filter ancillary related columns, e.g., the phase in simulations
        of periodic sources.
        """
        logger.info('Filtering %d event times according to the GTIs...', len(time_))
        mask = self.gti_mask(time_)
        time_ = time_[mask]
        logger.info('Done, %d entries remaining.', len(time_))
        return time_, mask

    def remove_bti(self, bti_start, bti_stop):
        """
        Update Good Time Intervals (GTIs) by removing Bad Time Intervals (BTIs)
        from them.
 
        Args:
            bti_start (numpy.ndarray): Array of start times for BTIs.
            bti_stop (numpy.ndarray): Array of stop times for BTIs.
 
        Returns:
            Tuple[numpy.ndarray, numpy.ndarray]: Arrays of start and stop times for
            the new GTIs.
        """
        gti_start, gti_stop = self.start_mets(), self.stop_mets()
        self.validate_intervals(bti_start, bti_stop)
        helper = xGTIListMergerHelper()
        
        def merge_rule():
            """ BTI are treated as negative GTIs
            """
            return helper.in_old_gti and not helper.in_new_gti
        
        return helper.combine_intervals(gti_start, gti_stop, bti_start,
                                        bti_stop, merge_rule)
    
    def merge_gti(self, new_gti_start, new_gti_stop):
        """
        Update Good Time Intervals (GTIs) by performing the *intersection* with 
        the given GTIs.
 
        Args:
            new_gti_start (numpy.ndarray): Array of start times for GTIs.
            new_gti_stop (numpy.ndarray): Array of stop times for GTIs.
 
        Returns:
            Tuple[numpy.ndarray, numpy.ndarray]: Arrays of start and stop times 
            for the intersection GTIs.
        """
        gti_start, gti_stop = self.start_mets(), self.stop_mets()
        self.validate_intervals(new_gti_start, new_gti_stop)
        helper = xGTIListMergerHelper()
        
        def merge_rule():
            return helper.in_old_gti and helper.in_new_gti
        
        return helper.combine_intervals(gti_start, gti_stop, new_gti_start,
                                        new_gti_stop, merge_rule)
 
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


class xGTIListMergerHelper:
    """ Helper class encapsulating most of the logic for combining two sets of  
    GTIs. Provides internal flags for keeping track of whether we are inside
    the OLD/NEW GTIs and the logic for switching state based on edge detection.
    """

    class EdgeType(IntEnum):
        """
        Enumeration representing the types of edges that can occur in Good
        Time Intervals (GTIs) and Bad Time Intervals (BTIs).

        Members:
            OLD_GTI_START: Beginning of a good time interval.
            OLD_GTI_STOP: End of a good time interval.
            NEW_GTI_START: Beginning of a bad time interval.
            NEW_GTI_STOP: End of a bad time interval.
        """
        OLD_GTI_START = 0
        OLD_GTI_STOP = 1
        NEW_GTI_START = 2
        NEW_GTI_STOP = 3


    EDGE_PRIORITY = {
        EdgeType.NEW_GTI_STOP:  0,
        EdgeType.OLD_GTI_STOP:  1,
        EdgeType.NEW_GTI_START: 2,
        EdgeType.OLD_GTI_START: 3,
    }

    def __init__(self):
        """
        """
        self.in_old_gti = False
        self.in_new_gti = False
    

    def _error_msg(self, edge_type):
        return  f'Unexpected edge type {edge_type.name} when '\
                f'"in_old_gti" is {self.in_old_gti} and '\
                f'"in_new_gti" is {self.in_new_gti}'

    def update_state(self, edge):
        """ Change the internal state based on the given edge type. Rise an 
        appropriate error in case of erroneous transitions (e.g. encoutering a 
        OLD_GTI_START when already inside a OLD_GTI). 
        """
        if edge == self.EdgeType.OLD_GTI_START:
            if self.in_old_gti:
                raise UnexpectedEdgeType(self._error_msg(edge))
            self.in_old_gti = True
        elif edge == self.EdgeType.OLD_GTI_STOP:
            if not self.in_old_gti:
                raise UnexpectedEdgeType(self._error_msg(edge))
            self.in_old_gti = False
        elif edge == self.EdgeType.NEW_GTI_START:
            if self.in_new_gti:
                raise UnexpectedEdgeType(self._error_msg(edge))
            self.in_new_gti = True
        elif edge == self.EdgeType.NEW_GTI_STOP:
            if not self.in_new_gti:
                raise UnexpectedEdgeType(self._error_msg(edge))
            self.in_new_gti = False
    
    def sweep_intervals(self, edge_times, edge_types, is_active):
        """
        Generic sweep-line engine for interval algebra.

        Intervals are half-open: [start, stop)

        The user provided function is_active defines when the resulting GTIs 
        should be considered active.

        Args
        ----------
        edge_times : ndarray
        edge_types : ndarray
        is_active : callable() -> bool
            Returns True when output GTI should be active.

        Returns
        -------
        start, stop : ndarray
        """

        # Sort edges by (time, priority)
        order = numpy.lexsort((
            [self.EDGE_PRIORITY[t] for t in edge_types],
            edge_times
        ))
        edge_times = edge_times[order]
        edge_types = edge_types[order]

        start = []
        stop = []

        prev_active = False

        for t, etype in zip(edge_times, edge_types):
            self.update_state(etype)
            active = is_active()

            if not prev_active and active:
                start.append(t)
            elif prev_active and not active:
                stop.append(t)

            prev_active = active

        if prev_active:
            raise RuntimeError('Unclosed interval at end of sweep')

        return numpy.array(start), numpy.array(stop)

    def combine_intervals(self, old_start, old_stop, new_start, new_stop, 
                          rule):
        """ Combine the given intervals based on the given rule. This function 
        mostly prepare the input arrays for sweep_intervals().
        """
        edge_times = numpy.hstack((old_start, old_stop, new_start, new_stop))
        edge_types = numpy.asarray(
            [self.EdgeType.OLD_GTI_START] * len(old_start) +
            [self.EdgeType.OLD_GTI_STOP]  * len(old_stop)  +
            [self.EdgeType.NEW_GTI_START] * len(new_start) +
            [self.EdgeType.NEW_GTI_STOP]  * len(new_stop)
        )
        return self.sweep_intervals(edge_times, edge_types, rule)


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

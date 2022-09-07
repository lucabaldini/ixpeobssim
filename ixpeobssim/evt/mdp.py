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

"""MDP-related facilities.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.utils.misc import pairwise
from ixpeobssim.utils.math_ import weighted_average
from ixpeobssim.utils.logging_ import logger

# pylint: disable=invalid-name


def mdp99(mu, signal_counts, bkg_counts=0.):
    """Return the MDP at the 99% confidence level.

    Note that the function returns numpy.inf if the number of signal events is
    zero.

    Arguments
    ---------
    mu : float or array
        The effective modulation factor (i.e., weighted over the count spectrum)

    signal_counts : float or array
        The number of signal events.

    bkg_counts : float or array
        The number of background events.
    """
    if numpy.count_nonzero(signal_counts) == 0:
        return numpy.inf
    return 4.292 / mu * numpy.sqrt(signal_counts + bkg_counts) / signal_counts



class xMDPRecord:

    """Small utility class to keep track of the ingredients going into a generic
    MDP calculation.
    """

    def __init__(self, emin, emax, mu, signal_counts, bkg_counts=0):
        """Constructor.
        """
        self.emin = emin
        self.emax = emax
        self.mu = mu
        self.signal_counts = signal_counts
        self.bkg_counts = bkg_counts
        self.mdp = mdp99(mu, signal_counts, bkg_counts)

    @staticmethod
    def empty(emin, emax):
        """Create an empty xMDPRecord object.

        This is useful when summing records in a loop as a zero-initializer.
        Note that the MDP corresponding to this object is infinite, which is
        harmless, as it effectively means it will be overwritten by any
        sensible class instance.
        """
        return xMDPRecord(emin, emax, 0., 0., 0.)

    def scale(self, scale_factor):
        """Scale the record by a simple multiplicative factor.

        Basic rules: the number of signal and background counts are scaled by
        the target amount, the effective modulation factor is unchanged, as the
        underlying spectrum is assumed to be un-modified, and the MDP is
        recalculated.
        """
        self.signal_counts *= scale_factor
        self.bkg_counts *= scale_factor
        self.mdp = mdp99(self.mu, self.signal_counts, self.bkg_counts)

    def __add__(self, other):
        """Overloaded operator to support addition.

        Note that we require MDP records to have the same energy bounds in
        order to be added,
        """
        assert self.emin == other.emin
        assert self.emax == other.emax
        signal_counts = self.signal_counts + other.signal_counts
        bkg_counts = self.bkg_counts + other.bkg_counts
        mu = weighted_average((self.mu, other.mu), (self.signal_counts, other.signal_counts))
        return self.__class__(self.emin, self.emax, mu, signal_counts, bkg_counts)

    def __str__(self):
        """String formatting.
        """
        text = '%.2f--%.2f keV: %10.1f' % (self.emin, self.emax, self.signal_counts)
        if self.bkg_counts > 0:
            text += '/%10.1f ' % self.bkg_counts
        text += ' counts, effective mu = %.3f, MDP = %.2f%%' % (self.mu, 100. * self.mdp)
        return text



class xMDPTable:

    """Small utility class to store a set MDP values evaluated in energy bins.

    This is essentially an ordered list of MDP records.

    TODO: allow calculating the sum row.
    """

    def __init__(self, observation_time):
        """Constructor.
        """
        self.observation_time = float(observation_time)
        self.source = None
        self.rows = []

    @staticmethod
    def empty(observation_time, ebinning):
        """Create an empty xMDPTable object.

        Like the corresponding function for the xMDPRecord class, this is useful
        when summing tables in a loop as a zero-initializer.

        Arguments
        ---------
        observation_time : float
            The duration of the observation for the MDP calculation

        ebinning : array_like
            The energy binning for the MDP table
        """
        table = xMDPTable(observation_time)
        for emin, emax in pairwise(ebinning):
            table.rows.append(xMDPRecord.empty(emin, emax))
        return table

    def set_source(self, source):
        """Attach a source to a given MPD table.

        This is implemented in this horrible way because a count spectrum is
        all you need to create an MDP table (i.e., you don't need an actual
        source), and in fact this is the way we operate in our
        applications---see xCountSpectrum.build_mdp_table(). But, after the
        fact, it might be useful to have access to the source parameters that
        the count spectrum corresponds to. We might want to review the
        implementation of the whole thing.
        """
        self.source = source

    def add_row(self, emin, emax, mu, signal_counts, bkg_counts=0):
        """Add a row to the MDP table.
        """
        row = xMDPRecord(emin, emax, mu, signal_counts, bkg_counts)
        self.rows.append(row)

    def num_rows(self):
        """Return the number of rows in the table.
        """
        return len(self.rows)

    def broadband_record(self):
        """Return a xMDPRecord with the broadband values.

        If the table has a single row, this is returning the row itself, with no
        additional calculations, otherwise the proper weighted average of the
        table rows is calculated.
        """
        assert self.num_rows() > 0
        if self.num_rows() == 1:
            return self.rows[-1]
        emin = self.rows[0].emin
        emax = self.rows[-1].emax
        _mu = [row.mu for row in self.rows]
        _sig = [row.signal_counts for row in self.rows]
        _bkg = [row.bkg_counts for row in self.rows]
        # Note that if there are no counts the effective modulation factor
        # is irrelevant and we set it to zero---this will happen if one tries
        # to print a table full of empty records, but other than that should be
        # relatively harmless.
        try:
            mu = weighted_average(_mu, _sig)
        except ZeroDivisionError:
            mu = 0.
        signal_counts = sum(_sig)
        bkg_counts = sum(_bkg)
        return xMDPRecord(emin, emax, mu, signal_counts, bkg_counts)

    def total_counts(self):
        """Return the total number of counts in the table.
        """
        return sum([row.signal_counts + row.bkg_counts for row in self.rows])

    def scale(self, scale_factor):
        """Scale the record by a simple multiplicative factor.

        See the corresponding xMDPRecord.scale() method for details.
        """
        for row in self.rows:
            row.scale(scale_factor)

    def __add__(self, other):
        """Overloaded operator for MDP table addition.
        """
        assert self.observation_time == other.observation_time
        assert self.num_rows() == other.num_rows()
        table = self.__class__(self.observation_time)
        for row1, row2 in zip(self.rows, other.rows):
            table.rows.append(row1 + row2)
        return table

    def __str__(self):
        """String formatting.
        """
        text = ''
        if self.source is not None:
            text += '%s\n' % self.source
        text += 'MDP table for %.1f s observation time\n' % self.observation_time
        for row in self.rows:
            text += '%s\n' % row
        if self.num_rows() > 1:
            text += '%s\n' % self.broadband_record()
        return text

    def save_ascii(self, file_path, **kwargs):
        """Save the table in ascii format.
        """
        logger.info('Writing MDP table to %s...', file_path)
        with open(file_path, 'w') as output_file:
            output_file.write('%s\n\n%s' % (kwargs, self))
        logger.info('Done.')

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

"""Miscellaneaous utilities.
"""

from __future__ import print_function, division

from itertools import tee

#pylint: disable=invalid-name


def pairwise(iterable):
    """Iterate over a binning vector.

    This will give you the bin edges for each bin, see
    https://stackoverflow.com/questions/5764782
    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def pairwise_enum(iterable):
    """Iterate over a binning vector keeping track of the indices.
    """
    return enumerate(pairwise(iterable))


def process_file_list(processing_function, file_list, *args, **kwargs):
    """Filter a series of files with a given processing function and return the
    list of the of the outputs from the processing, filtering out the None values.

    This is aimed at reducing the biolerplate code in the apps, where we typically
    make the same processing on a list of files and return the list of the
    processed files. Since for this application None is a signal for an existing
    file that does not get overwritten, filtering it out is helpful downstream.

    Args
    ----
    processing_function : callable
        The processing function to be called on the file---this should take the
        file path as a first argument.

    file_list : iterable
        The list of paths to the files to be processed.

    *args
        Positional arguments to be passed to the processing function.

    **kwargs
        Keyword arguments to be passed to the processing function.
    """
    output = []
    for file_path in file_list:
        ret = processing_function(file_path, *args, **kwargs)
        if ret is not None:
            output.append(ret)
    return output

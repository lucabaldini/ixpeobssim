#!/usr/bin/env python
#
# Copyright (C) 2021--2022, the ixpeobssim team.
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

"""xpevtstat app.
"""

from __future__ import print_function, division

import astropy
import numpy

from ixpeobssim.binning.base import xEventBinningBase
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.misc import pairwise

# pylint: disable=invalid-name

__description__ = \
"""Read a series of photon lists and print some basic statistics of the various
model components in the form of a simple text table.
"""


PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_ebinning(default_ebins=3)
PARSER.add_mc()



def xpevtstat(**kwargs):
    """Read the input file(s) and print the counts.
    """
    counts = {}
    ebinning = xEventBinningBase.make_energy_binning(**kwargs)
    for file_path in kwargs.get('filelist'):
        event_file = xEventFile(file_path)
        energy = event_file.energy_data(kwargs.get('mc'))
        roi_table = event_file.roi_table
        # If the input file has no ROI_TABLE, we only have one (unknown) source.
        if len(roi_table) == 0:
            n, _ = numpy.histogram(energy, ebinning)
            counts['Unknown'] = n
        # Otherwise, we loop over the ROI table.
        else:
            srcid = event_file.srcid_data()
            for source_id, source_name in roi_table.items():
                n, _ = numpy.histogram(energy[srcid == source_id], ebinning)
                if source_name in counts:
                    counts[source_name] += n
                else:
                    counts[source_name] = n
    total_counts = sum(counts.values())
    fractional_counts = {key: 100 * val / total_counts for key, val in counts.items()}
    col_names = ['Source name']
    col_names += ['%.2f--%.2f keV' % (emin, emax) for (emin, emax) in pairwise(ebinning)]
    col_types = [str] * len(ebinning)
    table = astropy.table.Table(names=col_names, dtype=col_types)
    for source_name in counts:
        values = ['%d (%.1f%%)' % (c, f) for c, f in zip(counts[source_name],
            fractional_counts[source_name])]
        table.add_row((source_name, *values))
    table.add_row(('Total', *['%d' % c for c in total_counts]))
    print('\n' * 3)
    print(table)



def main():
    """Main entry point.
    """
    xpevtstat(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
#
# Copyright (C) 2015, the ixpeobssim team.
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


__description__ = \
"""Tiny utility to query the CDS name resolver to attempt to retrieve
coordinate information for a given source.
"""


from astropy.coordinates import SkyCoord

from ixpeobssim.utils.logging_ import logger


def xpsrccoords(source_name):
    """
    """
    logger.info('Querying CDS name resolver for "%s"...' % source_name)
    coords = SkyCoord.from_name(source_name)
    print(coords.icrs)
    print(coords.galactic)
    logger.info('Done, bye!')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('source_name', type=str,
                        help='the name of the source')
    args = parser.parse_args()
    xpsrccoords(args.source_name)


if __name__ == '__main__':
    main()

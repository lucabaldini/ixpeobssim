#!/urs/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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

import os

from ixpeobssim import IXPEOBSSIM_DOC
from ixpeobssim.evt.event import xBinTableHDUMonteCarlo, xBinTableHDUGTI,\
    xBinTableHDURoiTable, xBinTableHDUTimeline, xBinTableHDUSpacecraftData,\
    xBinTableHDUOCTI, xLvl2PrimaryHDU, xBinTableHDUEvents
from ixpeobssim.evt.ixpesim import xBinTableHDUPhotons
from ixpeobssim.instrument.charging import xBinTableHDUCharging
from ixpeobssim.utils.logging_ import logger



def _title(text, char='='):
    """
    """
    return '%s\n%s\n\n' % (text, char * len(text))

def _subtitle(text):
    """
    """
    return _title(text, '-')

def _subsubtitle(text):
    """
    """
    return _title(text, '~')

def _format_extension(cls, name=None):
    """
    """
    if name is None:
        name = cls.NAME
    logger.info('Formatting extension %s', name)
    text = _subtitle('%s extension' % name)
    text += _subsubtitle('Header keywords')
    for key, value, comment in cls.HEADER_KEYWORDS:
        text += '* ``%s``: `%s`\n' % (key, comment)
    text += '\n'
    if not hasattr(cls, 'DATA_SPECS'):
        text += '\n'
        return text
    text += _subsubtitle('Columns')
    for key, fmt, units, comment in cls.DATA_SPECS:
        text += '* ``%s`` (%s): `%s` [%s]\n' % (key, fmt, comment, units)
    text += '\n'
    return text


def dump_format():
    """
    """
    logger.info('Dumping data format documentation...')
    file_path = os.path.join(IXPEOBSSIM_DOC, '_datafmt.rst')
    text = _format_extension(xLvl2PrimaryHDU, 'Primary')
    for cls in (xBinTableHDUEvents, xBinTableHDUGTI, xBinTableHDUMonteCarlo,
        xBinTableHDURoiTable, xBinTableHDUTimeline, xBinTableHDUSpacecraftData,
        xBinTableHDUOCTI, xBinTableHDUCharging):
        text += _format_extension(cls)
    logger.info('Writing output file %s...', file_path)
    with open(file_path, 'w') as output_file:
        output_file.write(text)
    logger.info('Done.')
    # 
    file_path = os.path.join(IXPEOBSSIM_DOC, '_phlistfmt.rst')
    text = _format_extension(xBinTableHDUPhotons)
    logger.info('Writing output file %s...', file_path)
    with open(file_path, 'w') as output_file:
        output_file.write(text)
    logger.info('Done.')



if __name__ == '__main__':
    dump_format()

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

from __future__ import print_function, division

import os
import shutil

from ixpeobssim import IXPEOBSSIM_DOC
from ixpeobssim.__version__ import TAG
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.os_ import cp


def package_docs():
    """
    """
    logger.info('Packaging the documentation for version %s...', TAG)
    build_dir = os.path.join(IXPEOBSSIM_DOC, '_build')
    src = os.path.join(build_dir, 'latex', 'ixpeobssim.pdf')
    dest = os.path.join(build_dir, 'ixpeobssim_%s.pdf' % TAG)
    cp(src, dest)
    src = os.path.join(build_dir, 'html')
    dest = os.path.join(build_dir, 'ixpeobssim_%s' % TAG)
    logger.info('Creating compressed html archive %s...', dest)
    shutil.make_archive(dest, 'zip', src)




if __name__ == '__main__':
    package_docs()

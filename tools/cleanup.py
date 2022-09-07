#!/usr/bin/env python
#
# * Copyright (C) 2015--2019, the ixpeobssim team.
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


import glob
import os
import shutil
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.system_ import cmd
from ixpeobssim import *


def cleanup(folder_path, patterns = ['*~', '*.pyc', '*.pyo'],
            cleanup_pycache=True):
    """Cleanup a folder.
    """
    logger.info('Cleaning up folder %s...' % folder_path)
    file_list = []
    for pattern in patterns:
        file_list += glob.glob(os.path.join(folder_path, pattern))
    for file_path in file_list:
        logger.info('Removing %s...' % file_path)
        os.remove(file_path)
    if cleanup_pycache:
        pycache_path = os.path.join(folder_path, '__pycache__')
        if os.path.exists(pycache_path):
            logger.info('Removing %s...' % pycache_path)
            shutil.rmtree(pycache_path)

def cleanup_dist():
    """Cleanup the distribution folder.
    """
    if os.path.exists(IXPEOBSSIM_DIST):
        logger.info('Removing %s altogether...' % IXPEOBSSIM_DIST)
        shutil.rmtree(IXPEOBSSIM_DIST)
    file_path = os.path.join(IXPEOBSSIM_ROOT, 'MANIFEST')
    if os.path.exists(file_path):
        logger.info('Removing %s...' % file_path)
        os.remove(file_path)

def cleanup_build():
    """
    """
    if os.path.exists(IXPEOBSSIM_BUILD):
        logger.info('Removing %s altogether...' % IXPEOBSSIM_BUILD)
        shutil.rmtree(IXPEOBSSIM_BUILD)

def cleanup_doc():
    """Cleanup the doc folder.
    """
    cmd('cd %s; make clean' % IXPEOBSSIM_DOC)

def cleanup_test_log():
    """Remove the unit test log file.
    """
    file_path = os.path.join(IXPEOBSSIM_TEST, 'unittest.log')
    if os.path.exists(file_path):
        logger.info('Removing %s...', file_path)
        os.remove(file_path)

def cleanup_setup():
    """
    """
    for folder in ('ixpeobssim.egg-info',):
        folder_path = os.path.join(IXPEOBSSIM_BASE, folder)
        if os.path.exists(folder_path):
            shutil.rmtree(folder_path)



if __name__ == '__main__':
    for folder_path in [
            IXPEOBSSIM_ROOT,
            IXPEOBSSIM_BENCHMARKS,
            IXPEOBSSIM_BIN,
            IXPEOBSSIM_BINNING,
            IXPEOBSSIM_CALDB,
            IXPEOBSSIM_CONFIG,
            IXPEOBSSIM_CORE,
            IXPEOBSSIM_DATA,
            IXPEOBSSIM_DIST,
            IXPEOBSSIM_DOC,
            IXPEOBSSIM_EVT,
            IXPEOBSSIM_EXAMPLES,
            IXPEOBSSIM_INSTRUMENT,
            IXPEOBSSIM_IRF,
            IXPEOBSSIM_IRFGEN,
            IXPEOBSSIM_NOTEBOOKS,
            IXPEOBSSIM_SANDBOX,
            IXPEOBSSIM_SRCMODEL,
            IXPEOBSSIM_TARGETS,
            IXPEOBSSIM_TEST,
            IXPEOBSSIM_TOOLS,
            IXPEOBSSIM_UTILS
    ]:
        cleanup(folder_path)
    cleanup_dist()
    cleanup_build()
    cleanup_doc()
    cleanup_test_log()
    cleanup_setup()

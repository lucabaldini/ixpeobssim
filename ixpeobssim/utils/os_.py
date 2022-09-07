#!/usr/bin/env python
# *********************************************************************
# * Copyright (C) 2015 Luca Baldini (luca.baldini@pi.infn.it)         *
# *                                                                   *
# * For the license terms see the file LICENSE, distributed           *
# * along with this software.                                         *
# *********************************************************************
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


"""Collection of os-related utilities.
"""


import os
import shutil

from ixpeobssim.utils.logging_ import logger, abort


def check_input_file(file_path, extension=None):
    """Make sure that an input file exists (and, optionally, has the right
    extension).

    Note that we abort the execution with no mercy if anything fails.
    """
    if not isinstance(file_path, str):
        abort('check_input_file() expecting a string, got %s' % file_path)
    if not os.path.exists(file_path):
        abort('Input file %s does not exists' % file_path)
    if not os.path.isfile(file_path):
        abort('Input file %s is not a file' % file_path)
    if extension is not None and not file_path.endswith('.%s' % extension):
        abort('Input file %s is not a .%s file' % (file_path, extension))


def check_output_file(file_path, suffix, overwrite=False, extension='fits'):
    """Small utility function to manage the I/O in the ixpeobssim applications.

    This basically build the path to the output file, given that to the input file,
    given a small set of rules. Among other things, the function verifies that
    the input file exists and has the right extension. If the output file
    exists and the overwite flag is not set, the function is returning None, as
    a mean for consumer functions to be aware that the aformentioned files should
    not be overwritten.

    Args
    ----
    file_path : str
        The path to the input file.

    suffix : str
        The suffix to be attached to the output file.

    overwrite : bool
        Flag to automatically overwite existing files.

    extension : str
        The extension of the input and output files.
    """
    assert extension is not None
    check_input_file(file_path, extension)
    output_file_path = file_path.replace('.%s' % extension, '_%s.%s' % (suffix, extension))
    if os.path.exists(output_file_path) and not overwrite:
        logger.info('Output file %s already exists.' % output_file_path)
        logger.info('Remove it or use "--overwrite True" to overwite it.')
        return None
    return output_file_path


def mkdir(dir_path):
    """Create a directory (unless it already exists).

    Return 0 upon succesfull operation, 1 otherwise.
    """
    if not os.path.exists(dir_path):
        logger.info('About to create folder %s...' % dir_path)
        try:
            os.makedirs(dir_path)
            logger.info('Folder succesfully created.')
            status = 0
        except Exception as e:
            logger.error('Could not create folder (%s)' % e)
            status = 1
        return status

def cp(source, dest, create_tree = False):
    """Copy a file.

    Return 0 upon succesfull operation, 1 otherwise.
    """
    logger.info('About to copy %s to %s...' % (source, dest))
    dest_folder = os.path.dirname(dest)
    if not os.path.exists(dest_folder) and create_tree:
        mkdir(dest_folder)
    try:
        if os.path.isdir(source):
            shutil.copytree(source, dest)
        else:
            shutil.copy(source, dest)
        logger.info('File succesfully copied.')
        status = 0
    except Exception as e:
        logger.error('Could not copy file (%s)' % e)
        status = 1
    return status

def mv(source, dest):
    """Move a file.

    Return 0 upon succesfull operation, 1 otherwise.
    """
    logger.info('About to move %s to %s...' % (source, dest))
    try:
        shutil.move(source, dest)
        logger.info('File succesfully copied.')
        status = 0
    except Exception as e:
        logger.error('Could not move file (%s)' % e)
        status = 1
    return status

def rm(file_path):
    """ Remove a file.

    Return 0 upon succesfull operation, 1 otherwise.
    """
    logger.info('About to remove file %s...' % file_path)
    if not os.path.exists(file_path):
        logger.info('File is not there, giving up...')
        return 0
    try:
        os.remove(file_path)
        logger.info('File succesfully removed.')
        status = 0
    except Exception as e:
        logger.error('Could not remove file (%s)' %  e)
        status = 1
    return status

def rmdir(dir_path):
    """ Remove an entire (empty or non empty) folder.
    """
    logger.info('About to remove folder %s...' % dir_path)
    try:
        shutil.rmtree(dir_path)
        logger.info('Folder succesfully removed.')
        status = 0
    except Exception as e:
        logger.error('Could not remove folder (%s)' %  e)
        status = 1
    return status

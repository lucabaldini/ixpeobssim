# Copyright (C) 2015--2019, the ixpeobssim team.
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

"""ixpeobssim: an X-ray polarimetry simulation framework.
"""

from __future__ import print_function, division

import os

import ixpeobssim.utils.os_
from ixpeobssim.utils.logging_ import logger

PACKAGE_NAME = 'ixpeobssim'

# Basic folder structure of the package.
IXPEOBSSIM_ROOT = os.path.abspath(os.path.dirname(__file__))
IXPEOBSSIM_BASE = os.path.abspath(os.path.join(IXPEOBSSIM_ROOT, os.pardir))
IXPEOBSSIM_BENCHMARKS = os.path.join(IXPEOBSSIM_ROOT, 'benchmarks')
IXPEOBSSIM_BIN = os.path.join(IXPEOBSSIM_ROOT, 'bin')
IXPEOBSSIM_BINNING = os.path.join(IXPEOBSSIM_ROOT, 'binning')
IXPEOBSSIM_BKG = os.path.join(IXPEOBSSIM_ROOT, 'bkg')
IXPEOBSSIM_BKG_DATA = os.path.join(IXPEOBSSIM_BKG, 'data')
IXPEOBSSIM_BUILD = os.path.join(IXPEOBSSIM_BASE, 'build')
IXPEOBSSIM_CALDB = os.path.join(IXPEOBSSIM_ROOT, 'caldb')
IXPEOBSSIM_CONFIG = os.path.join(IXPEOBSSIM_ROOT, 'config')
IXPEOBSSIM_CONFIG_ASCII = os.path.join(IXPEOBSSIM_CONFIG, 'ascii')
IXPEOBSSIM_CONFIG_FITS = os.path.join(IXPEOBSSIM_CONFIG, 'fits')
IXPEOBSSIM_CONFIG_REG = os.path.join(IXPEOBSSIM_CONFIG, 'reg')
IXPEOBSSIM_CORE = os.path.join(IXPEOBSSIM_ROOT, 'core')
IXPEOBSSIM_DIST = os.path.join(IXPEOBSSIM_BASE, 'dist')
IXPEOBSSIM_DOC = os.path.join(IXPEOBSSIM_BASE, 'docs')
IXPEOBSSIM_DOC_FIGURES = os.path.join(IXPEOBSSIM_DOC, 'figures')
IXPEOBSSIM_DOC_FIG_IRF = os.path.join(IXPEOBSSIM_DOC_FIGURES, 'IRF')
IXPEOBSSIM_DOC_FIG_MODELS = os.path.join(IXPEOBSSIM_DOC_FIGURES, 'models')
IXPEOBSSIM_DOC_FIG_OBSSIM = os.path.join(IXPEOBSSIM_DOC_FIGURES, 'obssim')
IXPEOBSSIM_DOC_FIG_MISC = os.path.join(IXPEOBSSIM_DOC_FIGURES, 'misc')
IXPEOBSSIM_DOC_FIG_TEST = os.path.join(IXPEOBSSIM_DOC_FIGURES, 'test')
IXPEOBSSIM_DOC_FIG_XSPEC = os.path.join(IXPEOBSSIM_DOC_FIGURES, 'xspec')
IXPEOBSSIM_EVT = os.path.join(IXPEOBSSIM_ROOT, 'evt')
IXPEOBSSIM_EXAMPLES = os.path.join(IXPEOBSSIM_ROOT, 'examples')
IXPEOBSSIM_INSTRUMENT = os.path.join(IXPEOBSSIM_ROOT, 'instrument')
IXPEOBSSIM_INSTRUMENT_DATA = os.path.join(IXPEOBSSIM_INSTRUMENT, 'data')
IXPEOBSSIM_IRF = os.path.join(IXPEOBSSIM_ROOT, 'irf')
IXPEOBSSIM_IRFGEN = os.path.join(IXPEOBSSIM_ROOT, 'irfgen')
IXPEOBSSIM_IRFGEN_DATA = os.path.join(IXPEOBSSIM_IRFGEN, 'data')
IXPEOBSSIM_NOTEBOOKS = os.path.join(IXPEOBSSIM_BASE, 'notebooks')
IXPEOBSSIM_OBSDATA = os.path.join(IXPEOBSSIM_ROOT, 'obsdata')
IXPEOBSSIM_SANDBOX = os.path.join(IXPEOBSSIM_ROOT, 'sandbox')
IXPEOBSSIM_SRCMODEL = os.path.join(IXPEOBSSIM_ROOT, 'srcmodel')
IXPEOBSSIM_TARGETS = os.path.join(IXPEOBSSIM_ROOT, 'targets')
IXPEOBSSIM_TARGETS_DATA = os.path.join(IXPEOBSSIM_TARGETS, 'data')
IXPEOBSSIM_TEST = os.path.join(IXPEOBSSIM_BASE, 'tests')
IXPEOBSSIM_TEST_DATA = os.path.join(IXPEOBSSIM_TEST, 'data')
IXPEOBSSIM_TOOLS = os.path.join(IXPEOBSSIM_BASE, 'tools')
IXPEOBSSIM_UTILS = os.path.join(IXPEOBSSIM_ROOT, 'utils')
IXPEOBSSIM_XSPEC = os.path.join(IXPEOBSSIM_ROOT, 'xspec')

# This is the output directory.
#
# Note that we create the folder if it does not exist.
try:
    IXPEOBSSIM_DATA = os.environ['IXPEOBSSIM_DATA']
except KeyError:
    IXPEOBSSIM_DATA = os.path.join(os.path.expanduser('~'), 'ixpeobssimdata')
ixpeobssim.utils.os_.mkdir(IXPEOBSSIM_DATA)
IXPEOBSSIM_DATA_BENCHMARKS = os.path.join(IXPEOBSSIM_DATA, 'benchmarks')
ixpeobssim.utils.os_.mkdir(IXPEOBSSIM_DATA_BENCHMARKS)


# Folder and infrastructure for auxiliary files.
try:
    IXPEOBSSIM_AUXFILES = os.environ['IXPEOBSSIM_AUXFILES']
except KeyError:
    IXPEOBSSIM_AUXFILES = os.path.join(os.path.expanduser('~'), 'ixpeobssimauxfiles')
ixpeobssim.utils.os_.mkdir(IXPEOBSSIM_AUXFILES)


AUXFILES_REPO_URL = 'https://bitbucket.org/ixpesw/ixpeobssim_auxfiles/downloads/'


def auxfile_path(file_name):
    """Return the full (local) path to an auxiliary file.
    """
    return os.path.join(IXPEOBSSIM_AUXFILES, file_name)

def auxfile_url(file_name):
    """Return the full (remote) path to an auxiliary file.
    """
    return os.path.join(AUXFILES_REPO_URL, file_name)

def auxfiles_missing(*file_list):
    """Check whether a set of auxiliary files exists in the local tree.
    """
    missing_files = []
    for file_name in file_list:
        if not os.path.exists(auxfile_path(file_name)):
            missing_files.append(file_name)
    if len(missing_files) > 0:
        logger.warning('You are missing the following auxiliary file(s).')
        for file_name in missing_files:
            logger.warning(auxfile_url(file_name))
        logger.warning('Please download the files manually to %s', IXPEOBSSIM_AUXFILES)
        return True
    return False



# Version information.
IXPEOBSSIM_VERSION_FILE_PATH = os.path.join(IXPEOBSSIM_ROOT, '__version__.py')

def version_info():
    """Read the tag and build date straight from the appropriate file.

    Use this when you don't want to import the module (i.e., at release time,
    when the file is changed), so that you don't have to bother with
    reloading stuff.
    """
    with open(IXPEOBSSIM_VERSION_FILE_PATH, 'r', encoding='utf8') as input_file:
        tag = input_file.readline().split('=')[-1].strip(' \'\n')
        build_date = input_file.readline().split('=')[-1].strip(' \'\n')
    return tag, build_date


# Release notes file.
IXPEOBSSIM_RELEASE_NOTES_PATH = os.path.join(IXPEOBSSIM_DOC, 'release_notes.rst')

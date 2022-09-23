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

"""Versioning of various third-party packages.
"""

import sys

import numpy
import scipy
import astropy
import matplotlib
import regions
import skyfield

from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.packaging_ import retrieve_version, xPackageVersion


# pylint: disable=unused-import


# Since the Python wrapper to xspec is not always trivial to set up, we put
# this guard, here, to signal downstream whether pyxspec is indeed installed
# or not.
#
# You can check whether the Python bindings for xspec are installed by
#
# >>> from ixpeobssim.utils.environment import PYXSPEC_INSTALLED
# >>> if PYXSPEC_INSTALLED:
# >>>     import ixpeobssim.evt.xspec_ as xspec_
#
# If PYXSPEC_INSTALLED is false you should refrain from touching anything into
#ixpeobssim.evt.xspec_, or event try to import anything from in there.
try:
    import xspec
    PYXSPEC_INSTALLED = True
except ImportError:
    logger.warning('PyXSPEC is not installed, you will no be able to use it.')
    PYXSPEC_INSTALLED = False

# Retrieve the Python version.
PYTHON_VERSION = xPackageVersion(sys.version_info.major, sys.version_info.minor,
    sys.version_info.micro)

# Retrieve the version numbers for some of the most important third-party packages.
NUMPY_VERSION = retrieve_version(numpy)
SCIPY_VERSION = retrieve_version(scipy)
ASTROPY_VERSION = retrieve_version(astropy)
MATPLOTLIB_VERSION = retrieve_version(matplotlib)
SKYFIELD_VERSION = retrieve_version(skyfield)
REGIONS_VERSION = retrieve_version(regions)

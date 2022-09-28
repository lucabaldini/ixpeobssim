#!/usr/bin/env python
#
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

from setuptools import setup, find_packages
import os
import sys
import glob

from ixpeobssim import PACKAGE_NAME
from ixpeobssim.__version__ import TAG


author = 'The ixpeobssim team'
description = 'Simulation and analysis framework for the Imaging X-ray Polarimetry Explorer'
license = 'GNU General Public License v3 or later'
packages = find_packages(exclude=('tests',))
url = 'https://github.com/lucabaldini/ixpeobssim/'
classifiers = [
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: '
    'GNU General Public License v3 or later (GPLv3+)',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: Implementation :: CPython',
    'Topic :: Scientific/Engineering :: Astronomy',
    'Development Status :: 5 - Production/Stable'
]
dependencies =[
    'numpy',
    'matplotlib',
    'astropy',
    'scipy',
    'regions',
    'skyfield'
]
entry_points = {}
entry_points['console_scripts'] = []
scripts = glob.glob(os.path.join(os.getcwd(), 'ixpeobssim', 'bin', '*.py'))

for scriptname in scripts:
    name = os.path.basename(scriptname.replace('.py', ''))
    func = 'ixpeobssim.bin.' + name + ':main'
    entry_points['console_scripts'].append('{0} = {1}'.format(name, func))

_KWARGS = dict(name=PACKAGE_NAME,
               version=TAG,
               author=author,
               description=description,
               long_description = open('README.md', 'rb').read().decode('utf-8'),
               long_description_content_type = 'text/markdown',
               license=license,
               packages=packages,
               include_package_data=True,
               url=url,
               classifiers=classifiers,
               scripts=[],
               install_requires=dependencies,
               entry_points=entry_points)

setup(**_KWARGS)

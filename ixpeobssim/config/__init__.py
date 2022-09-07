#!/usr/bin/env python
#
# Copyright (C) 2018, the ixpeobssim team.
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

from ixpeobssim import IXPEOBSSIM_CONFIG
from ixpeobssim.utils.matplotlib_ import plt, save_all_figures
from ixpeobssim.utils.logging_ import abort


def config_file_path(model_name):
    """Return the full path to the configuration file for a given model.
    """
    return os.path.join(IXPEOBSSIM_CONFIG, '%s.py' % model_name)


def file_path_to_model_name(file_path):
    """Transform a file path (typically inside the ixpeobssim/config folder)
    into a model name, i.e., strip the .py extension from the file name.
    """
    # Horrible hack for Python 2.
    if file_path.endswith('.pyc'):
        file_path = file_path.rstrip('c')
    if not file_path.endswith('.py'):
        abort('%s not a path to a Python file' % file_path)
    return os.path.basename(file_path).replace('.py', '')


def display_source_model(options, function=None, *args, **kwargs):
    """Small convenience function to minimize boilerplate code.

    This takes as an imput the arguments from out custom source model option
    parser and a disply function and will do all the magic.

    Arguments
    ---------
    options : argparse Namespace
        The arguments from the source model custom option parser

    function : function
        A reference to the Python function plotting the relevant stuff

    *args : 
        optional additional arguments to pass to the display function

    **kwargs :
        optional additional keyword arguments to pass to the display function
    """
    import __main__
    print(__main__.ROI_MODEL)
    if function is None:
        function = __main__.display
    function(*args, **kwargs)
    if options.save:
        save_all_figures(options.outfolder)
    if not options.batch:
        plt.show()


def bootstrap_display(function=None, *args, **kwargs):
    """Fire up the source model custom option parser and display the model.
    """
    from ixpeobssim.utils.argparse_ import xSourceModelArgumentParser
    parser = xSourceModelArgumentParser()
    display_source_model(parser.parse_args(), function, *args, **kwargs)

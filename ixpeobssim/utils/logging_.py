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


"""Logging utilities, building on top of the python logging module.
"""

import logging
import sys


class TerminalColors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def _color(text, color):
    """Process a piece of tect to be printed out in color.
    """
    return '%s%s%s' % (color, text, TerminalColors.ENDC)

def _red(text):
    """Process a piece of text to be printed out in red.
    """
    return _color(text, TerminalColors.RED)

def _yellow(text):
    """Process a piece of text to be printed out in yellow.
    """
    return _color(text, TerminalColors.YELLOW)

def _green(text):
    """Process a piece of text to be printed out in green.
    """
    return _color(text, TerminalColors.GREEN)


logger = logging.getLogger('ixpeobssim')
logger.setLevel(logging.DEBUG)


class xTerminalFormatter(logging.Formatter):

    """Logging terminal formatter class.
    """

    def format(self, record):
        """Overloaded format method.
        """
        text = ('>>> %s' % record.msg)
        if len(record.args) > 0:
            text = text % record.args
        if record.levelno >= logging.ERROR:
            text = _red(text)
        elif record.levelno == logging.WARNING:
            text = _yellow(text)
        return text


""" Configure the main terminal logger.
"""
consoleHandler = logging.StreamHandler()
consoleHandler.setLevel(logging.DEBUG)
consoleHandler.setFormatter(xTerminalFormatter())
logger.addHandler(consoleHandler)


def abort(message=''):
    """Abort the execution (via a sys.exit) with a message.

    Use this with care, and opt for custom exceptions whenever possible.
    """
    logger.error(message)
    sys.exit('Cannot continue, abort.')


def startmsg():
    """Print the start message.
    """
    from ixpeobssim.__version__ import TAG, BUILD_DATE
    print('\n    Welcome to ixpeobssim %s (built on %s).\n' %\
              (TAG, BUILD_DATE))
    print('    Copyright (C) 2015--2023, the ixpeobssim team.\n\n    ixpeobssim comes with ABSOLUTELY NO WARRANTY.\n    This is free software, and you are welcome to redistribute it under certain\n    conditions. See the LICENSE file for details.\n\n    Visit https://bitbucket.org/ixpesw/ixpeobssim for more information.\n')



if __name__ == '__main__':
    startmsg()

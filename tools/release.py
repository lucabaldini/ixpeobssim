#!/usr/bin/env python
#
# * Copyright (C) 2015, the ixpeobssim team.
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


import time
import os

from ixpeobssim.utils.os_ import rm
from ixpeobssim.utils.system_ import cmd
from ixpeobssim.utils.logging_ import logger
from ixpeobssim import IXPEOBSSIM_VERSION_FILE_PATH, version_info,\
    IXPEOBSSIM_RELEASE_NOTES_PATH, IXPEOBSSIM_DIST, IXPEOBSSIM_ROOT


BUILD_DATE = time.strftime('%a, %d %b %Y %H:%M:%S %z')
TAG_MODES = ['major', 'minor', 'patch']


def updateVersionInfo(mode, dry_run=False):
    """ Update the __tag__.py module with the new tag and build date.
    """
    prevTag, prevBuildDate = version_info()
    logger.info('Previous tag was %s...' % prevTag)
    version, release, patch = [int(item) for item in prevTag.split('.')]
    if mode == 'major':
        version += 1
        release = 0
        patch = 0
    elif mode == 'minor':
        release += 1
        patch = 0
    elif mode == 'patch':
        patch += 1
    else:
        abort('Unknown release mode %s.' % mode)
    nextTag = '%s.%s.%s' % (version, release, patch)
    logger.info('Writing new tag (%s) to %s...' %\
                (nextTag, IXPEOBSSIM_VERSION_FILE_PATH))
    if not dry_run:
        outputFile = open(IXPEOBSSIM_VERSION_FILE_PATH, 'w')
        outputFile.writelines('TAG = \'%s\'\n' % nextTag)
        outputFile.writelines('BUILD_DATE = \'%s\'\n' % BUILD_DATE)
        outputFile.close()
    logger.info('Done.')
    return nextTag

def updateReleaseNotes(tag, dry_run=False):
    """ Write the new tag and build date on top of the release notes
    (which must be kept up to date during the release process).
    """
    title = '.. _release_notes:\n\nRelease notes\n=============\n\n'
    logger.info('Reading in %s...' % IXPEOBSSIM_RELEASE_NOTES_PATH)
    notes = open(IXPEOBSSIM_RELEASE_NOTES_PATH).read().strip('\n').strip(title)
    logger.info('Writing out %s...' % IXPEOBSSIM_RELEASE_NOTES_PATH)
    if not dry_run:
        outputFile = open(IXPEOBSSIM_RELEASE_NOTES_PATH, 'w')
        outputFile.writelines(title)
        outputFile.writelines('\n*ixpeobssim (%s) - %s*\n\n' % (tag, BUILD_DATE))
        outputFile.writelines(notes)
        outputFile.close()
    logger.info('Done.')

def tagPackage(mode, dry_run=False):
    """ Tag the package.

    This means:
    (*) hg pull/update to make sure we're not missing remote modification;
    (*) figure out the target tag and update the release.notes;
    (*) commit the modifications, tag and push.
    """
    cmd('git pull', verbose=True, dry_run=dry_run)
    cmd('git status', verbose=True, dry_run=dry_run)
    tag = updateVersionInfo(mode, dry_run)
    updateReleaseNotes(tag, dry_run)
    msg = 'Prepare for tag %s.' % tag
    cmd('git commit -a -m "%s"' % msg, verbose=True, dry_run=dry_run)
    cmd('git push', verbose=True, dry_run=dry_run)
    msg = 'tagging version %s' % tag
    cmd('git tag -a %s -m "%s"' % (tag, msg), verbose=True, dry_run=dry_run)
    cmd('git push --tags', verbose = True, dry_run=dry_run)
    cmd('git status', verbose = True, dry_run=dry_run)

def distsrc():
    """ Create a plain source distribution.
    """
    tag, buildDate = version_info()
    logger.info('Creating plain source distribution...')
    distDir = os.path.join(IXPEOBSSIM_DIST, 'src')
    srcLogFilePath = 'src.log'
    # Create the distribution.
    cmd('python setup.py sdist --dist-dir=%s --prune' % distDir,
        verbose=False, logFilePath=srcLogFilePath)
    # Cleanup.
    rm(srcLogFilePath)
    rm(os.path.join(IXPEOBSSIM_ROOT, 'MANIFEST'))
    logger.info('Done.')



if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-t', dest = 'tagmode', type = str, default = None,
                      help = 'The release tag mode %s.' % TAG_MODES)
    parser.add_option('-n', action = 'store_true', dest = 'dryrun',
                      help = 'Dry run (i.e. do not actually do anything).')
    parser.add_option('-s', action = 'store_true', dest = 'src',
                      help = 'Create a source distribution.')
    (opts, args) = parser.parse_args()
    if not opts.tagmode and not (opts.src):
        parser.print_help()
        parser.error('Please specify at least one valid option.')
    tag = None
    if opts.tagmode is not None:
        if opts.tagmode not in TAG_MODES:
            parser.error('Invalid tag mode %s (allowed: %s)' %\
                             (opts.tagmode, TAG_MODES))
        tagPackage(opts.tagmode, opts.dryrun)
    if opts.src and not opts.dryrun:
        distsrc()

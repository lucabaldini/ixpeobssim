#!/usr/bin/env python
#
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

"""xpstokessmear app.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.core.hist import xGpdMap2d, xHistogram1d
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.instrument.gpd import GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_Y
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.utils.misc import process_file_list
from ixpeobssim.utils.os_ import check_output_file

#pylint: disable=no-member, invalid-name

__description__ = \
"""Smear the event-by-event Stokes parameters with a gaussian probability density
function with zero average.

This is supposed to be a cheap shortcut to test the effect of the correction for
the spurious modulation in realistic analysis without having to simulate the
effect to start with---we assume that the actual correction gets the average
right, and we are just introducing the statistical fluctuations of the
correction, here, in such a way that the output Stokes parameters are not
properly normalized, as will happen for real level-2 files.

The smearing is done over a fixed grid of pixels, with sigma values different in
the center of the detector and on the borders, and in a single energy band.
While this is defintely not capturing the entire richness of the phenomenology
of the spurious modulation, it is a sensible zero-order approximation to test
possible biases of the correction on actual science analysis.

For reference, the spurious modulation maps in the IXPE CALDB are mapped onto
300 x 300 bins on the detector surface, and across 6 energy bands. The typical
errors on Q and U are 0.025 in the inner 50 pixel (or 2.5 mm) radius circle,
and 0.10 in the outer portion of the detector.
"""

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_argument('--innersigma', type=float, default=0.025,
    help='the standard deviation for the Q and U gaussian smearing in the center')
PARSER.add_argument('--outersigma', type=float, default=0.1,
    help='the standard deviation for the Q and U gaussian smearing on the border')
PARSER.add_argument('--innerradius', type=float, default=2.5,
    help='the inner radius in mm, defining the region of the high-statistics calibration')
PARSER.add_argument('--nside', type=int, default=300,
    help='the size of the binning grid in DETX and DETY')
PARSER.add_suffix('stokessmear')
PARSER.add_seed()
PARSER.add_boolean('--interactive', default=False,
    help='plot some diagnostic plot')
PARSER.add_overwrite()


def _create_maps(**kwargs):
    """Create the correction map.

    Note this underwent some refactoring in response to issue
    https://github.com/lucabaldini/ixpeobssim/issues/668
    Particularly, since now the GPD maps are no more square, and they are binned
    over the physical area of the GPD, rather than the old fiducial square,
    the inner and outer masks are now calculated directly in mm rather than
    in pixels.

    Since this is a functionality that is rarely used in practice, we might
    want to watch out for possible issue with the changes.
    """
    if kwargs.get('seed') is not None:
        numpy.random.seed(kwargs.get('seed'))
    logger.info('Creating the Q and U correction maps...')
    nside = kwargs.get('nside')
    inner_sigma = kwargs.get('innersigma')
    outer_sigma = kwargs.get('outersigma')
    # Create the mask for the inner region.
    y, x = numpy.ogrid[:nside, :nside]
    x0 = y0 = nside / 2. - 0.5
    # Calculate the conversion factor between pixels and mm, so that we can
    # calculate the inner mask.
    xscale = 2. * GPD_PHYSICAL_HALF_SIDE_X / nside
    yscale = 2. * GPD_PHYSICAL_HALF_SIDE_Y / nside
    # Note the cast to float is a terrible workaround for issue
    # https://github.com/lucabaldini/ixpeobssim/issues/608
    # triggered by numpy 1.22.0
    x = xscale * x.astype(float)
    y = yscale * y.astype(float)
    x0 *= xscale
    y0 *= yscale
    # Calculate the actual masks.
    inner_mask = (x - x0)**2. + (y - y0)**2. <= kwargs.get('innerradius')**2.
    outer_mask = numpy.logical_not(inner_mask)
    ninner = inner_mask.sum()
    nouter = outer_mask.sum()
    # Create the actual smearing maps.
    qmap = xGpdMap2d(nside , 'Q smearing sigma')
    qmap.content[inner_mask] = numpy.random.normal(0., inner_sigma, ninner)
    qmap.content[outer_mask] = numpy.random.normal(0., outer_sigma, nouter)
    umap = xGpdMap2d(nside , 'U smearing sigma')
    umap.content[inner_mask] = numpy.random.normal(0., inner_sigma, ninner)
    umap.content[outer_mask] = numpy.random.normal(0., outer_sigma, nouter)
    return qmap, umap


def _plot_maps(qmap, umap):
    """Plotting function for debugging purposes.
    """
    plt.figure('Q smearing map')
    qmap.plot()
    plt.figure('U smearing map')
    umap.plot()


def _process_file(file_path, **kwargs):
    """Process a single file.
    """
    output_file_path = check_output_file(file_path, kwargs.get('suffix'), kwargs.get('overwrite'))
    if output_file_path is None:
        return None
    interactive = kwargs.get('interactive')
    # Create the smearing maps...
    qmap, umap = _create_maps(**kwargs)
    # Diagnostic plots...
    if interactive:
        plt.figure('Q smearing map')
        qmap.plot()
        plt.figure('U smearing map')
        umap.plot()
    # Open the input event file and cache the necessary columns.
    event_file = xEventFile(file_path)
    detx, dety = event_file.det_position_data()
    q, u = event_file.stokes_data()
    # Calculate the smearing deltas.
    logger.info('Calculating the smearing deltas...')
    dq = qmap.find_bin_value(detx, dety)
    du = umap.find_bin_value(detx, dety)
    logger.info('Delta Q: mean = %.3e, rms = %.3e', dq.mean(), dq.std(ddof=1))
    logger.info('Delta U: mean = %.3e, rms = %.3e', du.mean(), du.std(ddof=1))
    # Mofify the columns and write the output file.
    evt_data = event_file.hdu_list['EVENTS'].data
    evt_data['Q'] = q + dq
    evt_data['U'] = u + du
    # Remove the PHI column, if necessary.
    col_names = [col.name for col in event_file.hdu_list['EVENTS'].columns]
    if 'PHI' in col_names:
        event_file.remove_columns('EVENTS', 'PHI')
    event_file.write(output_file_path, overwrite=kwargs.get('overwrite'))
    return output_file_path


def xpstokessmear(**kwargs):
    """Application entry point.
    """
    return process_file_list(_process_file, kwargs.get('filelist'), **kwargs)


def main():
    """main() entry point.
    """
    xpstokessmear(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()

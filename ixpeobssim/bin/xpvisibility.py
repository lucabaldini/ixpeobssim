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

from __future__ import print_function, division, absolute_import


__description__ = \
"""Simple IXPE visibility tool.
"""

import numpy
import matplotlib.dates
from astropy.coordinates.name_resolve import NameResolveError
from astropy.coordinates import SkyCoord

from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.instrument.traj import xIXPETrajectory
from ixpeobssim.utils.time_ import string_to_met_utc, met_to_num,\
    years_to_seconds, met_to_string, seconds_to_days
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, metplot
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline


def _histogram(vals, binning, weights=None):
    """Small convienience function to bin data.
    """
    if weights is not None:
        weights = weights.astype(numpy.double)
    hist, _ = numpy.histogram(vals, binning, weights=weights)
    return hist


def xpvisibility(**kwargs):
    """.
    """
    # Calculate the start and stop MET for the observation.
    start_met = string_to_met_utc(kwargs.get('startdate'), lazy=True)
    if kwargs.get('stopdate') is not None:
        stop_met = string_to_met_utc(kwargs.get('stopdate'), lazy=True)
    else:
        stop_met = start_met + years_to_seconds(1.)
    # Parse the source coordinates, trying to resolve the source name, if passed.
    srcname = kwargs.get('srcname')
    if srcname is not None:
        logger.info('Trying to resolve source name "%s"...', srcname)
        try:
            coords = SkyCoord.from_name(srcname).icrs
        except NameResolveError as e:
            abort(e)
        logger.info(coords)
        ra, dec = coords.ra.deg, coords.dec.deg
    else:
        ra, dec = kwargs['srccoords']

    # Hard-coded parameters controlling the internals of the calculation.
    # Mind there are two natural scales for binning the fine-grid data, namely
    # the period of the orbit and the average period of the SAA ephocs.
    # Depending on which one we pick, the SAA and Earth occultation lines will
    # look smooth---but since the actual GTIs have all sort of resonances
    # between the two there is no easy way to get that smooth, too.
    num_samples = 100000
    orbit_period = 5805.14
    average_saa_period = 6212.78
    days = seconds_to_days(stop_met - start_met)
    scale = 0.2
    bin_width = average_saa_period * max(int(scale * days), 1)

    # Calculate the relevant masks along the trajectory.
    trajectory = xIXPETrajectory()
    met = numpy.linspace(start_met, stop_met, num_samples)
    logger.info('Sampling the orbit in %d MET values...', num_samples)
    saa = trajectory.in_saa(met)
    occult = trajectory.target_occulted(met, ra, dec)
    good = numpy.logical_not(saa + occult)

    # Bin the data.
    binning = numpy.arange(start_met, stop_met, bin_width)
    num_bins = len(binning)
    logger.info('Averaging values over %d (%.2f ks wide) bins...', num_bins, bin_width / 1000.)
    hmet = _histogram(met, binning)
    hsaa = _histogram(met, binning, saa) / hmet
    hoccult = _histogram(met, binning, occult) / hmet
    hgood = _histogram(met, binning, good) / hmet
    good_frac = hgood.mean()
    visibility_frac = 1. - hoccult.mean()

    # Print out the results.
    logger.info('Target coordinates: R. A. %.6f, Dec. %.6f', ra, dec)
    logger.info('Start MET for the observation: %.3f s', start_met)
    logger.info('Stop MET for the observation: %.3f s', stop_met)
    duration = 0.001 * (stop_met - start_met)
    logger.info('Wall-clock duration of the observation: %.3f ks', duration)
    logger.info('Total GTI: %.3f ks (%.1f %%)', good_frac * duration, 100 * good_frac)
    logger.info('Average visibility efficiency: %.3f%%', 100. * visibility_frac)

    # Create the splines for the plots.
    metc = 0.5 * (binning[:-1] + binning[1:])
    k = 3
    ssaa = xInterpolatedUnivariateSpline(metc, hsaa, k=k)
    soccult = xInterpolatedUnivariateSpline(metc, hoccult, k=k)
    sgood = xInterpolatedUnivariateSpline(metc, hgood, k=k)

    # Plot things.
    plt.figure('IXPE visibility tool')
    title = 'Estimated visibility for [%.4f$^\\circ$, %.4f$^\\circ$]' % (ra, dec)
    if srcname is not None:
        title = '%s (%s)' % (title, srcname)
    plt.title(title)
    met_grid = numpy.linspace(start_met, stop_met, 250)
    metplot(met_grid, ssaa(met_grid), label='In SAA')
    metplot(met_grid, soccult(met_grid), label='Target occulted by Earth')
    metplot(met_grid, sgood(met_grid), label='On source (GTIs)')
    setup_gca(ymin=0., ymax=1.1, ylabel='Fractional time in a given mode', grids=True)
    plt.legend(loc='upper left')

    if kwargs.get('sun'):
        sun_color = 'red'
        ax_left = plt.gca().twinx()
        setup_gca(ymin=0., ymax=200., ylabel='Sun pitch angle [$^\\circ$]')
        ax_left.spines['right'].set_color(sun_color)
        ax_left.tick_params(axis='y', colors=sun_color)
        ax_left.yaxis.label.set_color(sun_color)
        thin_line_fmt = dict(color=sun_color, ls='dashed', lw=0.75)
        label_fmt = dict(color=sun_color, rotation=90., va='center', size='small')
        met = numpy.linspace(start_met, stop_met, 250)
        sun_angle = trajectory.target_sun_angle(met, ra, dec)
        ax_left.plot(met_to_num(met), sun_angle, color=sun_color)
        # Draw some cute small Suns along the pitch angle :-)
        delta = 0.05 * (stop_met - start_met)
        x = numpy.linspace(start_met + delta, stop_met - delta, 9)
        y = trajectory.target_sun_angle(x, ra, dec)
        plt.plot(met_to_num(x), y, '.', marker='$\\odot$', ms=9., color=sun_color,
            label='Sun pitch angle')
        plt.legend(loc='upper right')
        ax_left.axhline(trajectory.min_sun_angle, **thin_line_fmt)
        ax_left.axhline(trajectory.max_sun_angle, **thin_line_fmt)
        epochs = trajectory.sun_constraint_epochs(start_met, stop_met, ra, dec)
        _format = lambda met: met_to_string(met1).split('T')[0]
        for met1, met2 in epochs:
            t1 = met_to_num(met1)
            t2 = met_to_num(met2)
            plt.gca().axvspan(t1, t2, alpha=0.10, color='gray', hatch='//')
            #plt.axvline(t1, **thin_line_fmt)
            #plt.axvline(t2, **thin_line_fmt)
            #if met1 > start_met:
            #    plt.text(t1, 90., _format(met1), ha='left', **label_fmt)
            #if met2 < stop_met:
            #    plt.text(t2, 90., _format(met2), ha='right', **label_fmt)
        epochs = trajectory.complement_epochs(epochs, start_met, stop_met)
        logger.info('Viewing periods (due to Sun constraints):')
        for i, (met1, met2) in enumerate(epochs):
            logger.info('[%02d] %s--%s', i + 1, met_to_string(met1), met_to_string(met2))

    plt.tight_layout()



"""Command-line switches.
"""
from ixpeobssim.utils.argparse_ import xArgumentParser

PARSER = xArgumentParser(description=__description__)
PARSER.add_target_source()
PARSER.add_startdate()
PARSER.add_stopdate()
PARSER.add_boolean('--sun', default=True, help='plot the Sun constraints')
PARSER.add_overwrite()


def main():
    kwargs = PARSER.parse_args().__dict__
    xpvisibility(**kwargs)



if __name__ == '__main__':
    main()
    plt.show()

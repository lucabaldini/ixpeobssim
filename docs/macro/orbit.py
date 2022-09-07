#!/urs/bin/env python
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

from __future__ import print_function, division

import sys

import numpy
import matplotlib.dates

from ixpeobssim.instrument.traj import xIXPETrajectory, xSAABoundary, xTLE
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, DEFAULT_COLORS
from ixpeobssim.utils.matplotlib_ import save_all_figures
from ixpeobssim.utils.time_ import string_to_met_utc, met_to_num
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.utils.logging_ import logger
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim import IXPEOBSSIM_DOC_FIG_MISC

if sys.flags.interactive:
    plt.ion()



def plot(start_date = '2021-06-01', duration=86400, save=False):
    """
    """
    start_met = string_to_met_utc(start_date, lazy=True)
    stop_met = start_met + duration
    met = numpy.linspace(start_met, stop_met, 5000)
    trajectory = xIXPETrajectory()
    saa = xSAABoundary()
    geo_fmt = dict(xlabel='Longitude [deg]', ylabel='Latitude [deg]',
                   xmin=-180., xmax=180., grids=True)

    print('%s\n%s' % xTLE.lines())
    altitude = 600.
    mean_motion = xTLE._mean_motion(altitude)
    period = 86400. / mean_motion
    print('Mean motion = %.5f, period = %.5f' % (mean_motion, period))

    # Calculate the trajectory and altitude.
    lon, lat = trajectory.position(met)
    alt = trajectory.elevation(met)
    # Split the trajectory into segements to avoid annoying horizontal lines
    # in the plot when passing from 180 to -180 in longitude.
    idx = numpy.where(numpy.diff(lon) < 0.)[0] + 1
    lon = numpy.split(lon, idx)
    lat = numpy.split(lat, idx)

    # Plain satellite orbit.
    plt.figure('IXPE trajectory')
    for _lon, _lat in zip(lon, lat):
        plt.plot(_lon, _lat, color=DEFAULT_COLORS[0])
    plt.plot([lon[0][0], lon[-1][-1]], [lat[0][0], lat[-1][-1]], 'o')
    setup_gca(**geo_fmt)

    # Altitude
    plt.figure('IXPE altitude')
    plt.plot(met_to_num(met), alt)
    setup_gca(ylabel='Altitude [km]', ymin=599.6, ymax=600.4, grids=True)
    locator = matplotlib.dates.AutoDateLocator()
    formatter = matplotlib.dates.ConciseDateFormatter(locator)
    plt.gca().xaxis.set_major_formatter(formatter)

    # SAA
    plt.figure('IXPE SAA polygon')
    saa.plot()
    setup_gca(ymin=-35., **geo_fmt)

    # Satellite orbit with SAA.
    plt.figure('IXPE trajectory SAA')
    for _lon, _lat in zip(lon, lat):
        plt.plot(_lon, _lat, color=DEFAULT_COLORS[0])
        mask = saa.contains(_lon, _lat)
        plt.plot(_lon[mask], _lat[mask], color='lightgray')
    saa.plot()
    setup_gca(ymin=-31., ymax=6., **geo_fmt)
    axins = plt.gca().inset_axes([0.575, 0.35, 0.4, 0.4])
    axins.set_xlim(-78., 0.)
    axins.set_ylim(-0.55, 0.55)
    axins.set_xticklabels('')
    axins.set_yticklabels('')
    for _lon, _lat in zip(lon, lat):
        axins.plot(_lon, _lat, color=DEFAULT_COLORS[0])
        mask = saa.contains(_lon, _lat)
        axins.plot(_lon[mask], _lat[mask], color='lightgray')
    fmt = dict(facecolor='orange', edgecolor='black', alpha=0.5)
    patch = matplotlib.patches.PathPatch(saa, **fmt)
    axins.add_patch(patch)
    plt.plot(*numpy.hsplit(saa.vertices, 2), 'o', color=fmt.get('edgecolor'))

    plt.gca().indicate_inset_zoom(axins)
    axins.grid(which='both')

    # SAA passages
    epochs = trajectory.saa_epochs(start_met, start_met + 1000000, 700)
    delta = numpy.array([t2 - t1 for (t1, t2) in epochs])
    min_ = delta.min()
    max_ = delta.max()
    mean = delta.mean()
    rms = delta.std(ddof=1)
    plt.figure('IXPE SAA epochs')
    hist = xHistogram1d(numpy.linspace(780., 840., 75)).fill(delta)
    hist.plot()
    setup_gca(ymax=17., grids=True, xlabel='SAA epoch duration [s]', ylabel='Entries per bin')

    pipeline.xpvisibility(srcname='Crab', startdate='2021-01-01')

    if save:
        save_all_figures(IXPEOBSSIM_DOC_FIG_MISC, ('pdf', 'png'))



if __name__ == '__main__':
    plot(save=True)

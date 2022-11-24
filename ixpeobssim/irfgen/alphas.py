#!/usr/bin/env python
#
# Copyright (C) 2020--2021, the ixpeobssim team.
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

"""Small module providing simulation tools for the alpha particles emitted by
the Be window.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.irfgen.gpd import WINDOW_BE_THICKNESS, WINDOW_AL_THICKNESS
from ixpeobssim.irfgen.gpd import ABSORPTION_GAP_THICKNESS
from ixpeobssim.instrument.gpd import GPD_PHYSICAL_HALF_SIDE_X
from ixpeobssim.core.geometry import xPoint, xLine, xRay
from ixpeobssim.irfgen.astar import load_alpha_stopping_power_data
from ixpeobssim.irfgen.astar import load_dme_alpha_stopping_power_data
from ixpeobssim.irfgen.constants import BE_DENSITY, AL_DENSITY
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.math_ import modulo_2pi


WINDOW_HALF_SIDE = 4. #cm
ASIC_HALF_SIDE = GPD_PHYSICAL_HALF_SIDE_X / 10. #cm
# Setup the geometry---all dimensions in cm.
Z_ASIC = 0.
Z_AL_BOT = Z_ASIC + ABSORPTION_GAP_THICKNESS
Z_BE_BOT = Z_AL_BOT + WINDOW_AL_THICKNESS
Z_BE_TOP = Z_BE_BOT + WINDOW_BE_THICKNESS
logger.info('Basic detector geometry...')
logger.info('ASIC at %.6f cm', Z_ASIC)
logger.info('Window Al bottom at %.8f cm', Z_AL_BOT)
logger.info('Window Be bottom at %.8f cm', Z_BE_BOT)
logger.info('Window Be top at %.8f cm', Z_BE_TOP)

BE_ASTAR_DATA = load_alpha_stopping_power_data('Be', BE_DENSITY)
AL_ASTAR_DATA = load_alpha_stopping_power_data('Al', AL_DENSITY)
DME_ASTAR_DATA = load_dme_alpha_stopping_power_data()



def _inside(val, minimum, maximum):
    """Convenience function returning True if a value is within a min and max.
    """
    return val >= minimum and val <= maximum

def _zinside(val):
    """
    """
    return _inside(val, Z_ASIC, Z_AL_BOT)

def _xyinside(val):
    """
    """
    return _inside(val, -ASIC_HALF_SIDE, ASIC_HALF_SIDE)

def gas_cell_intersection_points(ray):
    """Calculate the intersecting points of a given ray with the gas cell.

    This is a horrible list of conditions, and should be
    """
    w = ASIC_HALF_SIDE
    east = ray.yintersect(ASIC_HALF_SIDE)
    west = ray.yintersect(-ASIC_HALF_SIDE)
    north = ray.xintersect(ASIC_HALF_SIDE)
    south = ray.xintersect(-ASIC_HALF_SIDE)
    top = ray.zintersect(Z_AL_BOT)
    bottom = ray.zintersect(Z_ASIC)
    points = []
    for point in [east, west]:
        if _zinside(point.z) and _xyinside(point.x):
            points.append(point)
    for point in [north, south]:
        if _zinside(point.z) and _xyinside(point.y):
            points.append(point)
    for point in [top, bottom]:
        if _xyinside(point.x) and _xyinside(point.y):
            points.append(point)
    num_points = len(points)
    # If there are no valid intersection points, return None and handle this
    # downstream.
    if num_points == 0:
        return None
    # At this pooint the intersection points must be exactly 2.
    assert num_points == 2
    # Since the alpha is travelling downward, make sure that the z of the first
    # point is higher than the z of the second one.
    p1, p2 = points
    if p1.z < p2.z:
        points = [p2, p1]
    return points


def generate_decays(alpha_energy=5., num_alphas=1000000, plot=True, energy_cut=0.1):
    """Generate the initial alpha particles in the Be window volume.

    This is returning a list of xRay object and a numpy array of energies
    describing the alpha particles that emerge from the window in a direction
    that is intersecting the fiducial cylinder enclosing the GPD active gas
    volume.

    This first part of the program is vectorized, as only a small fraction of
    the particles survive the initial cut, and the whole thing would be
    excruciatingly slow, if stuffed into a plain Python loop.
    """
    logger.info('Generating %d initial alpha decays...', num_alphas)
    # Generate the z decay poition within the Be window.
    z = numpy.random.uniform(Z_BE_BOT, Z_BE_TOP, size=num_alphas)
    # Generate the z cosine director uniformly in the half-sphere.
    zdir = numpy.random.random(size=num_alphas)
    # Calculate the path length in the Berillium
    be_path_length = (z - Z_BE_BOT) / zdir
    logger.info('Average decay z: %.8f cm', z.mean())
    # Calculate the range of an alpha particle in the Be.
    be_range = BE_ASTAR_DATA.csda_range(alpha_energy)
    logger.info('Alpha range in Be: %.8f um', be_range * 1.e4)

    # Initial debug plots.
    if plot:
        plt.figure('Conversion point')
        binning = numpy.linspace(Z_BE_BOT, Z_BE_TOP, 100)
        h = xHistogram1d(binning, xlabel='z conversion point [cm]').fill(z)
        h.plot()
        plt.figure('Path length')
        binning = numpy.linspace(0., 0.025, 100)
        h = xHistogram1d(binning, xlabel='Path length in Be [cm]').fill(be_path_length)
        h.plot()
        plt.vlines(be_range, 0., num_alphas, ls='dashed', color='gray')

    # Throw away all the events that would not make it out of the Be window.
    mask = be_path_length < be_range
    n = mask.sum()
    # Calculate the fraction of alpha escaping the Be Window
    # (note the factor of 0.5 is due to the fact that we are generating over
    # 2pi steradians.)
    frac = 100. * 0.5 * n / num_alphas
    logger.info('%d alpha(s) escaping the Be (%.2f%%)...', n, frac)
    # Select the particles that survive.
    be_path_length = be_path_length[mask]
    z = z[mask]
    zdir = zdir[mask]
    # Calculate the energy of the surviving alpha particles at the Be exit.
    be_energy_profile = BE_ASTAR_DATA._energy_profile(alpha_energy)
    energy = be_energy_profile(be_path_length)

    if plot:
        plt.figure('Be energy profile')
        be_energy_profile.plot(grids=True)

    # Now the energy loss in the Al layer.
    al_path_length = WINDOW_AL_THICKNESS / zdir
    al_energy_loss = AL_ASTAR_DATA.stopping_power(energy) * al_path_length
    logger.info('Average energy loss in the Al layer: %.2e MeV', al_energy_loss.mean())
    mask = energy > al_energy_loss
    n = mask.sum()
    frac = 100. * 0.5 * n / num_alphas
    logger.info('%d alpha(s) escaping the Al (%.2f%%)...', n, frac)
    z = z[mask]
    zdir = zdir[mask]
    energy = energy[mask]

    # Apply the minimum energy cut.
    mask = energy > energy_cut
    n = mask.sum()
    frac = 100. * 0.5 * n / num_alphas
    logger.info('%d alpha(s) above %.2f MeV (%.2f%%)...', n, energy_cut, frac)
    z = z[mask]
    zdir = zdir[mask]
    energy = energy[mask]

    # More debug plots.
    if plot:
        plt.figure('Conversion point (gas volume)')
        binning = numpy.linspace(Z_BE_BOT, Z_BE_TOP, 100)
        h = xHistogram1d(binning, xlabel='z conversion point [cm]').fill(z)
        h.plot()
        plt.figure('Alpha energy (gas volume)')
        binning = numpy.linspace(0., alpha_energy, 100)
        h = xHistogram1d(binning, xlabel='Alpha energy [MeV]').fill(energy)
        h.plot()

    # Generate the starting decay position within the Be window.
    x = numpy.random.uniform(-WINDOW_HALF_SIDE, WINDOW_HALF_SIDE, size=n)
    y = numpy.random.uniform(-WINDOW_HALF_SIDE, WINDOW_HALF_SIDE, size=n)
    phi = 2. * numpy.pi * numpy.random.random(size=n)

    # Additional geometrical cut to throw away all the rays that do not intercept
    # the circle enclosing the active volume of the gas.
    # See https://math.stackexchange.com/questions/543496 for the basic math.
    r = ASIC_HALF_SIDE * numpy.sqrt(2.)
    d = numpy.sqrt(x**2. + y**2)
    rho = r / d
    b = d * rho * numpy.sqrt(1. - rho**2.)
    phi0 = numpy.pi + numpy.arctan2(y, x)
    delta_phi = numpy.arccos(b / r)
    mask = numpy.logical_or((rho > 1.), (abs(modulo_2pi(phi - phi0)) < delta_phi))
    n = mask.sum()
    frac = 100. * 0.5 * n / num_alphas
    logger.info('%d alpha(s) within fiducial cylinder (%.2f%%)...', n, frac)

    # Prepare the output rays.
    x = x[mask]
    y = y[mask]
    z = numpy.full(x.shape, Z_AL_BOT)
    theta = numpy.degrees(numpy.arccos(-zdir[mask]))
    phi = numpy.degrees(phi[mask])
    energy = energy[mask]

    # Even more debug plots.
    if plot:
        plt.figure('Alpha energy (fiducial cylinder)')
        binning = numpy.linspace(0., alpha_energy, 100)
        h = xHistogram1d(binning, xlabel='Alpha energy [MeV]').fill(energy)
        h.plot()

    if plot:
        plt.show()

    rays = []
    for _x, _y, _z, _theta, _phi in zip(x, y, z, theta, phi):
        p = xPoint(_x, _y, _z)
        ray = xRay(p, _theta, _phi)
        rays.append(ray)
    return tuple(rays), energy



def plot_rectangle(width, height, center=(0., 0.), **kwargs):
    """Facility to plot a generic rectangle.
    """
    x0, y0 = center
    x1 = x0 - 0.5 * width
    x2 = x0 + 0.5 * width
    y1 = y0 - 0.5 * height
    y2 = y0 + 0.5 * height
    plt.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], **kwargs)



def plot_event(ray, enery, track_start=None, track_end=None):
    """Rudimentary single-event display.
    """
    p0 = ray.origin
    p1 = ray.zintersect(Z_ASIC)
    ray_fmt = dict(ls='dashed', color='orange', lw=1.)
    window_fmt = dict(ls='solid', color='gray', lw=1.)
    gas_fmt = dict(ls='solid', color='blue', lw=1.)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, num='Event display', figsize=(9., 9.))

    plt.sca(ax3)
    plot_rectangle(2. * WINDOW_HALF_SIDE, 2. * WINDOW_HALF_SIDE, **window_fmt)
    plot_rectangle(2. * ASIC_HALF_SIDE, 2. * ASIC_HALF_SIDE, **gas_fmt)
    plt.axis([-5., 5., -5., 5.])
    plt.gca().set_aspect('equal')
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    plt.plot([p0.x, p1.x], [p0.y, p1.y], **ray_fmt)
    for point in [track_start, track_end]:
        plt.plot([track_start.x, track_end.x], [track_start.y, track_end.y], 'o')

    plt.sca(ax1)
    plot_rectangle(2. * WINDOW_HALF_SIDE, ABSORPTION_GAP_THICKNESS,
                   center=(0, 0.5 * ABSORPTION_GAP_THICKNESS), **window_fmt)
    plot_rectangle(2. * ASIC_HALF_SIDE, ABSORPTION_GAP_THICKNESS,
                   center=(0, 0.5 * ABSORPTION_GAP_THICKNESS), **gas_fmt)
    plt.axis([-5., 5., -5., 5.])
    plt.gca().set_aspect('equal')
    plt.xlabel('x [cm]')
    plt.ylabel('z [cm]')
    plt.plot([p0.x, p1.x], [p0.z, p1.z], **ray_fmt)
    for point in [track_start, track_end]:
        plt.plot([track_start.x, track_end.x], [track_start.z, track_end.z], 'o')

    plt.sca(ax4)
    plot_rectangle(2. * WINDOW_HALF_SIDE, ABSORPTION_GAP_THICKNESS,
                   center=(0, 0.5 * ABSORPTION_GAP_THICKNESS), **window_fmt)
    plot_rectangle(2. * ASIC_HALF_SIDE, ABSORPTION_GAP_THICKNESS,
                   center=(0, 0.5 * ABSORPTION_GAP_THICKNESS), **gas_fmt)
    plt.axis([-5., 5., -5., 5.])
    plt.gca().set_aspect('equal')
    plt.xlabel('y [cm]')
    plt.ylabel('z [cm]')
    plt.plot([p0.y, p1.y], [p0.z, p1.z], **ray_fmt)
    for point in [track_start, track_end]:
        plt.plot([track_start.y, track_end.y], [track_start.z, track_end.z], 'o')



class AlphaTrack:

    def __init__(self, start, end, deposited_energy):
        """Constructor.
        """
        self.start = start
        self.end = end
        self.deposited_energy = deposited_energy
        self.length = self.start.distance_to(self.end)

    def projected_length(self):
        """Return the length of the track, projected on the xy plane.
        """
        dx = self.end.x - self.start.x
        dy = self.end.y - self.start.y
        return numpy.sqrt(dx**2. + dy**2.)

    def roi_size(self, padding=0.025, pixel_area=2.17e-5):
        """
        """
        dx = abs(self.end.x - self.start.x) + 2. * padding
        dy = abs(self.end.y - self.start.y) + 2. * padding
        return int(dx * dy / pixel_area)

    def track_size(self, width=0.05, pixel_area=2.17e-5):
        """
        """
        return int(self.length * width / pixel_area)



def simulate(alpha_energy=5.0, num_alphas=100000, plot=False, energy_cut=0.1):
    """Run the actual simulation.
    """
    track_list = []

    # Setup the energy loss for the DME.
    dme_energy_loss = DME_ASTAR_DATA.energy_loss(alpha_energy)

    # Create the list of particle direction and energy for the alpha particles
    # that emerge from the window with a direction intersecting the fiducial
    # cylinder enclosing the active gas cell.
    logger.info('Propagating tracks...')
    for ray, energy in zip(*generate_decays(alpha_energy, num_alphas, plot)):
        points = gas_cell_intersection_points(ray)
        # If the event does not hit the active gas volume, there is nothing to do.
        if points is None:
            continue
        # Unpack the enter and exit pointo into/from the gas active volume.
        track_start, track_end = points
        # Calculate the energy at the entry point in the active gas volume.
        energy -= dme_energy_loss(energy, ray.origin.distance_to(track_start))
        # Does the particle die before reaching the active gas volume?
        if energy <= energy_cut:
            continue
        # Build the actual track and energy loss.
        track_length = track_start.distance_to(track_end)
        deposited_energy = dme_energy_loss(energy, track_length)
        # If the particle looses all the remaining energy, it dies within the
        # active volume, and we have to go back and re-calculate the track
        # endpoint.
        if energy - deposited_energy < 1e-6:
            track_length = DME_ASTAR_DATA.csda_range(energy)
            deposited_energy = energy
            track_end = track_start.move(track_length, ray)

        track_list.append(AlphaTrack(track_start, track_end, deposited_energy))

        if plot:
            plot_event(ray, energy, track_start, track_end)
            plt.show()

    n = len(track_list)
    frac = 100. * 0.5 * n / num_alphas
    logger.info('Done, %d good tracks generated (%.2f%%).', n, frac)
    return track_list


def generate_distributions(alpha_energy=5.0, num_alphas=10000000):
    """
    """
    track_list = simulate(alpha_energy, num_alphas)
    deposited_energy = [track.deposited_energy for track in track_list]
    length = [track.length for track in track_list]
    projected_length = [track.projected_length() for track in track_list]
    roi_size = [track.roi_size() for track in track_list]
    track_size = [track.track_size() for track in track_list]
    binning = numpy.linspace(0., alpha_energy, 100)
    hde = xHistogram1d(binning, xlabel='Deposited energy [MeV]').fill(deposited_energy)
    binning = numpy.linspace(0., 2., 100)
    htl = xHistogram1d(binning, xlabel='Track length [cm]').fill(length)
    binning = numpy.linspace(0., 2., 100)
    hptl = xHistogram1d(binning, xlabel='Projected track length [cm]').fill(projected_length)
    binning = numpy.linspace(0., 40000, 100)
    hrs = xHistogram1d(binning, xlabel='ROI size [pixels]').fill(roi_size)
    binning = numpy.linspace(0., 10000, 100)
    hts = xHistogram1d(binning, xlabel='Track size [pixels]').fill(track_size)
    return hde, htl, hptl, hrs, hts


def plot_distributions(alpha_energy=5.0, num_alphas=1000000):
    """
    """
    hde, htl, hptl, hrs, hts = generate_distributions(alpha_energy, num_alphas)
    plt.figure('Deposited energy')
    hde.plot()
    plt.figure('Track length')
    htl.plot()
    plt.figure('Projected track length')
    hptl.plot()
    plt.figure('ROI size')
    hrs.plot()
    plt.figure('Track size')
    hts.plot()



def plot_raindrops(alpha_energy=5.0, num_alphas=2000, max_roi_size=5000):
    """
    """
    plot_rectangle(2. * ASIC_HALF_SIDE, 2. * ASIC_HALF_SIDE)
    plt.axis([-1., 1., -1., 1.])
    plt.gca().set_aspect('equal')
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    for track in simulate(alpha_energy, num_alphas):
        if track.roi_size() > max_roi_size:
            fmt = dict(lw=1., ls='dashed', color='lightgray')
        else:
            fmt = dict(ls='solid', color='gray')
        plt.plot([track.start.x, track.end.x], [track.start.y, track.end.y], **fmt)



if __name__ == '__main__':
    plot_distributions()
    plt.figure('Event raindrops')
    plot_raindrops()
    plt.show()

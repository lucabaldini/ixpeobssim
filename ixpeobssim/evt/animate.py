#!/usr/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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


"""Animation utilities.
"""

from __future__ import print_function, division

import os
import numbers

import numpy
import matplotlib.patches
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.image import AxesImage

from ixpeobssim import IXPEOBSSIM_CONFIG_FITS
from ixpeobssim.core.fitsio import xFITSImageBase
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.utils.astro import angular_separation
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.utils.units_ import degrees_to_arcsec, arcmin_to_degrees



class xMovingCircle:

    """Small convenience class representing a time-dependent circular patch.

    The moving patch is initialized given a set of (t, x, y) points, referred
    to the center of the target image. The interpolation is done via
    interpolating splines of order 1.

    Arguments
    ---------
    points : iterable of three-elements (t, x, y) tuples
        The points defining the path of the circle as a function of time.
        (The time is in seconds and the xy position is expressed as the delta
        with respect to the center of the image in arcminutes.)

    radius : the circle radius in arcseconds
    """

    def __init__(self, points, radius=25., kxy=1, kr=1):
        """Constructor.
        """
        if isinstance(radius, numbers.Number):
            radius = numpy.full(len(points), radius)
        else:
            radius = numpy.array(radius)
        t = numpy.array([p[0] for p in points])
        assert radius.shape == t.shape
        self.tmin = t[0]
        self.tmax = t[-1]
        dx = numpy.array([p[1] for p in points])
        dy = numpy.array([p[2] for p in points])
        self._splinex = xInterpolatedUnivariateSpline(t, dx, k=kxy)
        self._spliney = xInterpolatedUnivariateSpline(t, dy, k=kxy)
        self._spliner = xInterpolatedUnivariateSpline(t, radius, k=kr)

    def __call__(self, t):
        """Overloaded method.

        This is returning the x and y positions, in arcmin from the image center,
        calculated on a given time grid.
        """
        return self._splinex(t), self._spliney(t), self._spliner(t)

    def frame_positions(self, interval=100):
        """Return a set of frame positions, referred to the center of the target
        image, for given animation interval.

        The positions are in the form (dx, dy, radius), where dx and dy are the
        distances to the image center along the two orthogonal axes, in arcmin,
        and radius is the circle radius in arcseconds---this has to match the
        signature of the xSkyAnimation.show_roi() method, which is used as the
        updating slot of the animation.

        See https://matplotlib.org/stable/api/animation_api.html for more
        information about the matplotlib animation APIs.

        Arguments
        ---------
        interval : the animation interval in ms.
        """
        # Calculate the time grid for the interpolation---this is essentially
        # a constant-step grid spanning the entire time interval defined by the
        # input point, with the step defined by the target animation interval.
        # (Note the conversion from ms to s.)
        t = numpy.arange(self.tmin, self.tmax, 1.e-3 * interval)
        dx, dy, r = self(t)
        return list(zip(dx, dy, r))

    def event_mask(self, t, ra, dec, ra0, dec0):
        """Return a boolean mask for a series of events.

        Arguments
        ---------
        t : array_like
            The array of event times.

        ra : array_like
            The array of the event R.A. positions.

        dec : array_like
            The array of the event Dec. positions.

        ra0 : float
            The R.A. position of the target image center.

        dec0 : float
            The Dec. position of the target image center.
        """
        dx, dy, r = self(t)
        ra0, dec0 = ra0 + arcmin_to_degrees(dx), dec0 + arcmin_to_degrees(dy)
        angsep = degrees_to_arcsec(angular_separation(ra, dec, ra0, dec0))
        return angsep <= r



class xSkyAnimation:

    """Small convenience class to describe an animation in sky coordinates.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        self.image = xFITSImageBase(file_path)
        self.ra0, self.dec0 = self.image.center()

    @staticmethod
    def _find_last_axes_image():
        """Find the last children AxesImage object in the current axes.

        See https://stackoverflow.com/questions/25505341
        for the origin of this horrible hack.
        """
        for child in plt.gca().get_children():
            if isinstance(child, AxesImage):
                return child

    def plot(self, **kwargs):
        """Plot the underlying image.

        Note that, in addition to dispatch the plot() call to the underlying
        xFITSImageBase object, this is overriding the return value,
        returning an AxesImage object, rather than the standard figure. This is
        needed to properly build the animation video.
        """
        kwargs.setdefault('stretch', 'log')
        self.image.plot(**kwargs)
        return self._find_last_axes_image()

    def add_roi(self, roi, **kwargs):
        """Add the ROI.
        """
        delta_ra, delta_dec, radius = roi
        kwargs.setdefault('color', 'black')
        kwargs.setdefault('alpha', 0.50)
        kwargs.setdefault('animated', True)
        ra = self.ra0 + arcmin_to_degrees(delta_ra)
        dec = self.dec0 +  arcmin_to_degrees(delta_dec)
        x, y = self.image.wcs.wcs_world2pix(ra, dec, 0)
        radius /= degrees_to_arcsec(abs(self.image.cdelt2()))
        try:
            patch = matplotlib.patches.Annulus((x, y), 1000., 1000. - radius, **kwargs)
        except AttributeError:
            patch = matplotlib.patches.Circle((x, y), radius, **kwargs)
        plt.gca().add_patch(patch)
        return patch

    def show_frame(self, roi=(0., 0., 25.), **kwargs):
        """Draw a circular ROI with a given position and radius.

        Arguments
        ---------
        roi : three-element tuple
            The position and radius of the circular ROI in the form of a
            three-element tuple (dx, dy, radius).
        """
        plt.clf()
        im = self.plot()
        patch = self.add_roi(roi, **kwargs)
        return im, patch

    def run(self, roi, interval=100, figsize=(12., 12.)):
        """Run the actual animation.

        See https://matplotlib.org/stable/api/animation_api.html for more
        information about the matplotlib animation APIs.

        See also https://holypython.com/how-to-save-matplotlib-animations-the-ultimate-guide/

        Arguments
        ---------
        interval : the animation interval in ms.
        """
        plt.figure(figsize=figsize)
        frames = roi.frame_positions(interval)
        anim = FuncAnimation(plt.gcf(), self.show_frame, frames=frames, blit=True,
            interval=interval, repeat=False, save_count=len(frames))
        return anim

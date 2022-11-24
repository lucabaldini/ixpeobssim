#!/urs/bin/env python
#
# Copyright (C) 2018--2022, the ixpeobssim team.
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


import numpy

from ixpeobssim.instrument.sc import dithering_pattern
from ixpeobssim.utils.units_ import arcmin_to_degrees, degrees_to_arcsec
from ixpeobssim.instrument.gpd import rotate_detxy
from ixpeobssim.utils.logging_ import logger

# The focal length has been updated form the nominal 4-m value to the measured
# value of 3997 mm based on issue #336.
FOCAL_LENGTH = 3997.


def fiducial_backscal(half_side_x, half_side_y):
    """Calculate the BACKSCALE value for a give fiducial rectangle in detector
    coordinates.

    This function was introduced in response to the issue
    https://github.com/lucabaldini/ixpeobssim/issues/668
    in place of the old FIDUCIAL_BACKSCAL constant.

    Arguments
    ---------
    half_side_x : float
        The half side of the fiducial rectangle along the x coordinate in mm.

    half_side_y : float
        The half side of the fiducial rectangle along the y coordinate in mm.
    """
    dx = numpy.degrees(2. * half_side_x / FOCAL_LENGTH)
    dy = numpy.degrees(2. * half_side_y / FOCAL_LENGTH)
    return degrees_to_arcsec(dx) * degrees_to_arcsec(dy)

def _sky_to_gpd_naive(ra, dec, ra_pnt, dec_pnt):
    """Convert an array of (ra, dec) positions in the sky to an array
    of (detx, dety) positions onto the GPD reference frame.

    As the name (and the prepending _) suggest, this is a simplified version
    of the transformation, not taking into account a number of effects,
    including the relative rotation of the DUs, the roll angle and the dithering
    of the observatory. It is still useful to provide an implementation that is
    agnostic about the specific identifier of the DU, the roll angle and the
    event time---which is necessary for calculting the dithering information.
    Note that this function is called twice in the actual sky_to_gpd().

    Arguments
    ---------
    ra : float or array
        The ra position in the sky in decimal degrees

    dec : float or array
        The dec position in the sky in decimal degrees

    ra_pnt : float or array
        The pointing ra in decimal degrees

    dec_pnt : float or array
        The pointing dec in decimal degrees
    """
    cos_dec = numpy.cos(numpy.radians(dec_pnt))
    detx = - FOCAL_LENGTH * numpy.radians(ra - ra_pnt) * cos_dec
    dety = FOCAL_LENGTH * numpy.radians(dec - dec_pnt)
    return detx, dety

def parse_dithering_kwargs(**kwargs):
    """Parse the keyword arguments related to the dithering.
    """
    if not kwargs.get('dithering'):
        return None
    return (kwargs.get('ditherampl'), kwargs.get('ditherpa'),
        kwargs.get('ditherpx'), kwargs.get('ditherpy'))

def _dithering_delta(time_, dither_params, degrees=True):
    """Return the angular offset of the observatory boresight wrt to the
    nominal pointing position at a given array of event times and for a given
    set of dithering parameters.

    Note that the ditheing pattern expresses the orientation of the observatory
    with respect to the nominal pointing direction, so that if you want
    to turn this into the actual displacement in the sky of a given photon in
    the event list due to the dithering you should take the opposite and take
    cosine of the declination into account downstream.

    By default the return value is expressed in decimal degrees.

    Arguments
    ---------
    time_ : float or array
        The time at which the angular offset should be calculated

    dither_params : tuple or list
        The complete set of dithering parameters, in the following order
        (ditherampl, ditherpa, ditherpx, ditherpy), in the same units they would
        be passed to xpobssim via command-line.

    degrees : bool
        If true (default) the agular offset in ra and dec is converted in
        decimal degrees.
    """
    logger.info('Calculating the dithering angular offset...')
    logger.info('A = %.3f arcmin, pa = %.3f s, px = %.3f s, py = %.3f s', *dither_params)
    dithering = dithering_pattern(*dither_params)
    # Compute the dithering offset based on time---note this is in arcmin!
    delta_ra, delta_dec = dithering(time_)
    if degrees:
        # Convert from arcmin to degrees.
        delta_ra = arcmin_to_degrees(delta_ra)
        delta_dec = arcmin_to_degrees(delta_dec)
    return delta_ra, delta_dec


def apply_dithering(time_, ra_pnt, dec_pnt, dither_params=None):
    """Apply dithering correction directly to pointing direction in celestial coordinates.
    In this way the correction gets propagated down the chain and the appropriate IRFs are used.
    See convolve_event_list in mma.py

    Arguments
    ---------

    time : float or array
        The event time

    ra_pnt : float or array
        The pointing ra in decimal degrees

    dec_pnt : float or array
        The pointing dec in decimal degrees

    dither_params : (amplitude, pa, px, py) tuple, optional
        The parameters for the dithering of the observatory. The dithering is
        not applied if this is set to None.
    """
    # Apply the dithering.
    if dither_params is None:
        return numpy.full(time_.shape, ra_pnt), numpy.full(time_.shape, dec_pnt)
    delta_ra, delta_dec = _dithering_delta(time_, dither_params)
    ra_pnt += delta_ra / numpy.cos(numpy.radians(dec_pnt))
    dec_pnt += delta_dec
    return ra_pnt, dec_pnt


def sky_to_gpd(ra, dec, time_, ra_pnt, dec_pnt, du_id, roll_angle, dither_params=None):
    """Convert an array of (ra, dec) positions in the sky to an array
    of (x, y) positions onto the gpd reference frame. This involves essentially
    three different steps: first we project the sky coordinates onto the focal
    plane reference frame, then we apply the dithering pattern (if enabled)
    based on the event times and finally we rotate the coordinates according to
    the DU id and the telescope roll angle.

    Warning
    -------
    Change x and y to detx and dety!

    Arguments
    ---------
    ra : float or array
        The ra position in the sky in decimal degrees

    dec : float or array
        The dec position in the sky in decimal degrees

    time : float or array
        The event time

    ra_pnt : float or array
        The pointing ra in decimal degrees

    dec_pnt : float or array
        The pointing dec in decimal degrees

    du_id : int
        The detector unit id

    roll_angle : float
        The telescope roll angle in decimal degrees

    dither_params : (amplitude, pa, px, py) tuple, optional
        The parameters for the dithering of the observatory. The dithering is
        not applied if this is set to None.
    """
    # Project sky position to the "naive" reference frame
    detx, dety = _sky_to_gpd_naive(ra, dec, ra_pnt, dec_pnt)
    if dither_params is not None:
        # Apply the dithering.
        delta_ra, delta_dec = _dithering_delta(time_, dither_params)
        # Project the dithering onto the focal plane by calling the naive
        # transformation centered on the pointing of the observatory.
        dx, dy = _sky_to_gpd_naive(delta_ra, delta_dec, 0., 0.)
        detx -= dx
        dety -= dy
    # Rotate the xy coordinates according to du id and roll angle
    return rotate_detxy(detx, dety, du_id, roll_angle)


def _gpd_to_sky_naive(detx, dety, ra_pnt, dec_pnt):
    """Convert an array of (detx, dety) positions in detector coordinates to an
    array of (ra, dec) positions in the sky. This is supposed to be the inverse
    of _sky_to_gpd_naive(), so that the two operations, applied in series, should
    preserve the original data.

    Arguments
    ---------
    detx : float or array
        The detx position in detector coordinates in mm

    dety : float or array
        The dety position in detector coordinates in mm

    ra_pnt : float or array
        The pointing ra in decimal degrees

    dec_pnt : float or array
        The pointing dec in decimal degrees
    """
    cos_dec = numpy.cos(numpy.radians(dec_pnt))
    ra = ra_pnt - numpy.degrees(detx / FOCAL_LENGTH / cos_dec)
    dec = dec_pnt + numpy.degrees(dety / FOCAL_LENGTH)
    return ra, dec


def gpd_to_sky(detx, dety, time_, ra_pnt, dec_pnt, du_id, roll_angle, dither_params=None):
    """Convert an array of (detx, dety) positions in detector coordinates to an
    array of (ra, dec) positions in the sky. This is supposed to be the inverse
    of sky_to_gpd(), so that the two operations, applied in series, should
    preserve the original data.

    Arguments
    ---------
    detx : float or array
        The detx position in detector coordinates in mm

    dety : float or array
        The dety position in detector coordinates in mm

    time : float or array
        The event time

    ra_pnt : float or array
        The pointing ra in decimal degrees

    dec_pnt : float or array
        The pointing dec in decimal degrees

    du_id : int
        The detector unit id

    roll_angle : float
        The telescope roll angle in decimal degrees

    dither_params : (amplitude, pa, px, py) tuple, optional
        The parameters for the dithering of the observatory. The dithering is
        not applied if this is set to None.
    """
    detx, dety = rotate_detxy(detx, dety, du_id, roll_angle, inverse=True)
    ra, dec = _gpd_to_sky_naive(detx, dety, ra_pnt, dec_pnt)
    if dither_params is not None:
        # Apply the dithering.
        delta_ra, delta_dec = _dithering_delta(time_, dither_params)
        ra += delta_ra / numpy.cos(numpy.radians(dec_pnt))
        dec += delta_dec
    return ra, dec

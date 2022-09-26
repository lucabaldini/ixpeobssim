# Copyright (C) 2019, the ixpeobssim team.
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

"""Astronomical utilities.
"""


from __future__ import print_function, division

import operator
import numbers
import functools

from astropy import units, wcs
from astropy.coordinates import angle_utilities, SkyCoord
import numpy
import regions

from ixpeobssim.utils.environment import REGIONS_VERSION
from ixpeobssim.utils.logging_ import logger


# pylint: disable=invalid-name, no-member

REGION_COMPOUND_DICT = {'or': operator.or_, 'and': operator.and_, 'xor': operator.xor}
REGION_COMPOUND_OPERATORS = tuple(REGION_COMPOUND_DICT.keys())


def read_ds9(file_path):
    """Read a ds9 file from disk.

    This is a horrible hack to support multiple versions of the regions package,
    see https://bitbucket.org/ixpesw/ixpeobssim/issues/589

    Arguments
    ---------
    file_path : str
        The path to the actual file.
    """
    if REGIONS_VERSION <= '0.5':
        return regions.read_ds9(file_path)
    else:
        return regions.Regions.read(file_path, format='ds9')


def region_compound(*region_list, compound_mode='or'):
    """Create a region compound from an arbitrary number of astropy regions.
    """
    # If only one region is passed, there is nothing to do.
    if len(region_list) == 1:
        return region_list[0]
    assert compound_mode in REGION_COMPOUND_OPERATORS
    return functools.reduce(REGION_COMPOUND_DICT[compound_mode], region_list[1:], region_list[0])


def ds9_region_filter_xy(x, y, *region_list, compound_mode="or", invert=False, wcs_=None):
    """Check which the x and y sky coordinates of a series of events are inside
    the proper logical combination of a list of regions.

    Args
    ----
    x : array-like
        The input array of x coordinates (projection of sky coordinates).

    y : array
        The input array of y coordinates (projection of sky coordinates).

    region_list :
        A series of astropy pixel regions used to filter the events. (Note you
        can pass a list of sky regions, provided you also pass the proper
        wcs object to turn that into a pixel region, using the wcs_ argument.)

    compound_mode: string
        The logic by which to combine regions ('or', 'xor', or 'and').

    invert: bool
        A boolean flag that provides the mask that filters the inverse of
        the region list provided

    wcs_ : wcs object, optional
        If passed the code will attempt to turn a sky region into a pixel
        region based on the wcs information.

    Returns
    -------
    array-like
        The mask array corresponding to the input events.
    """
    reg = region_compound(*region_list, compound_mode=compound_mode)
    if wcs_ is not None:
        reg = reg.to_pixel(wcs_)
    mask = reg.contains(regions.PixCoord(x, y))
    if invert:
        mask = numpy.invert(mask)
    return mask


def ds9_region_filter_sky(ra, dec, wcs_, *region_list, compound_mode="or", invert=False):
    """Check which the ra and dec sky coordinates of a series of events are inside
    the proper logical combination of a list of regions, given the associated wcs.

    Args
    ----
    ra : array-like
        The input array of ra coordinates (projection of sky coordinates).

    dec : array
        The input array of dec coordinates (projection of sky coordinates).

    wcs_ : astropy.wcs.WCS object
        The wcs_ object mapping sky coordinates into pixels.

    region_list :
        A series of astropy sky regions used to filter the events.

    compound_mode: string
        The logic by which to combine regions ('or', 'xor', or 'and').

    invert: bool
        A boolean flag that provides the mask that filters the inverse of
        the region list provided

    Returns
    -------
    array-like
        The mask array corresponding to the input events.
    """
    reg = region_compound(*region_list, compound_mode=compound_mode)
    mask = reg.contains(SkyCoord(ra * units.degree, dec * units.degree), wcs_)
    if invert:
        mask = numpy.invert(mask)
    return mask


def angular_separation(ra, dec, ra0, dec0):
    """Compute the angular separation (in decimal degrees) of an array of
    coordinates wrt to a reference position in the sky.

    This replaces the old angular_distance() method, and is not using the
    astropy.SkyCoord facilities, which proved to be excruciatingly slow.
    The implementation uses the Vincenty formula, see
    https://en.wikipedia.org/wiki/Great-circle_distance
    as implemented in astropy.angle_utilities. (The only thing we are really
    doing here is to convert degrees to radians to degrees.)

    Args
    ----
    ra : array-like
        The input array of right ascension coordinates.

    dec : array-like
        The input array of declination coordinates.

    ra0 : float
        The reference ra coordinate.

    dec0 : float
        The reference dec coordinate.
    """
    args = [numpy.radians(val) for val in (ra, dec, ra0, dec0)]
    return numpy.degrees(angle_utilities.angular_separation(*args))


def square_sky_grid(nside, center, half_size):
    """Create a square, regular grid in ra and dec.

    Arguments
    ---------
    nside : int
        The number of points for either side of the grid

    center : (float, floar)
        The ra, dec center of the grid in decimal degrees.

    half_size : float
        The half-size of the grid in decimal degrees.
    """
    ra0, dec0 = center
    grid = numpy.linspace(-half_size, half_size, nside)
    ra = grid / numpy.cos(numpy.radians(dec0)) + ra0
    dec = grid + dec0
    return numpy.meshgrid(ra, dec)


def _calculate_crpix(nside):
    """Calculate the value of the TCRPX keyword for a given WCS size---we factor
    this small functionality in a separate function so that we can reuse it
    in wcs_columns_kwargs() and in build_wcs().

    We set TCRPX, corresponding to the input ra and dec, in the middle of the
    image, i.e., at (nside + 1) / 2. (e.g., if nside is 3, CRPIX is 2, which
    is in the middle of the [1, 2, 3] sequence, starting from 1 in the FITS
    convention.)
    """
    return 0.5 * (nside + 1)


def xy_columns_kwargs(ra0, dec0, nside, pixel_size, projection='TAN'):
    """
    Return a list of keyword arguments to be used for the X and Y columns (i.e.,
    those with a WCS associated.)

    This was added to fix some astropy warnings
    https://bitbucket.org/ixpesw/ixpeobssim/issues/523/
    and it definitely feels like the right thing to do. See the astropy docs
    for more details:
    https://docs.astropy.org/en/stable/io/fits/api/tables.html#astropy.io.fits.Column

    Arguments
    ---------
    ra0 : float
        The right ascension of the center of the wcs in decimal degrees.

    dec0 : float
        The declination of the center of the wcs in decimal degrees.

    nside : int
        The number of pixels in either dimension.

    pixel_size : float
        The WCS pixel size in decimal degrees.
    """
    crpx = _calculate_crpix(nside)
    xkwargs = dict(coord_type='RA---%s' % projection, coord_unit='deg',
        coord_ref_point=crpx, coord_ref_value=ra0, coord_inc=-pixel_size)
    ykwargs = dict(coord_type='DEC--%s' % projection, coord_unit='deg',
        coord_ref_point=crpx, coord_ref_value=dec0, coord_inc=pixel_size)
    return xkwargs, ykwargs


def set_xy_header_limits(hdu, nside):
    """Set the TLMIN and TLMAX keywords for the X and Y columns for a given hdu.

    I am somewhat surprised that this needs to be done by hand, but it is
    relevant, e.g., if one wants to properly display an event file with ds9.
    Here we are essentially looping over all the columns in the hdu and
    tweaking by hamd those that are literally called X and Y.
    """
    logger.info('Updating X and Y columns limits (TLMINn and TLMAXn)...')
    for i, col in enumerate(hdu.columns):
        if col.name in ('X', 'Y'):
            n = (i + 1)
            for key, val in (('TLMIN%d' % n, 1), ('TLMAX%d' % n, nside)):
                logger.debug('Setting %s to %d', key, val)
                hdu.header.set(key, val)


def build_wcs(ra0, dec0, nside, pixel_size, projection='TAN'):
    """Build a (square) WCS object programmatically.

    For reference, here are the basic WCS fields that we need to define,
    from https://idlastro.gsfc.nasa.gov/ftp/pro/astrom/aaareadme.txt

    * NAXIS: 2 element long vector giving dimensions of the images
    * CD: 2 x 2 array containing the astrometry parameters CD1_1 CD1_2
      in DEGREES/PIXEL
    * CDELT: 2 element vector giving physical increment at the reference pixel
    * CRPIX: 2 element vector giving X and Y coordinates of reference pixel
      (def = NAXIS/2) in FITS convention (first pixel is 1,1)
    * CRVAL: 2 element double precision vector giving R.A. and DEC of
      reference pixel in DEGREES
    * CTYPE: 2 element string vector giving projection types, default
      ['RA---TAN','DEC--TAN']
    * LONGPOLE: scalar longitude of north pole (default = 180)
    * LATPOLE: scalar giving native latitude of the celestial pole default=0)
    * PV2: Vector of parameters (PV2_1, PV2_2...) needed in some projections
    * DISTORT: Optional substructure giving distortion parameters. Currently
      only implemented for the Spitzer simple imaging polynomial (SIP) see
      http://fits.gsfc.nasa.gov/registry/sip.html

    Arguments
    ---------
    ra0 : float
        The right ascension of the center of the wcs in decimal degrees.

    dec0 : float
        The declination of the center of the wcs in decimal degrees.

    nside : int
        The number of pixels in either dimension.

    pixel_size : float
        The WCS pixel size in decimal degrees.
    """
    logger.debug('Building wcs object...')
    wcs_side = nside * pixel_size
    logger.debug('%d pixel(s) @ %.5f deg (%.3f deg image size)', nside, pixel_size, wcs_side)
    # The number of axes must be set from the start.
    wcs_ = wcs.WCS(naxis=2)
    # 2-element vector giving physical increment at the reference pixel.
    wcs_.wcs.cdelt = -pixel_size, pixel_size
    # 2-element vector giving X and Y coordinates of reference pixel.
    _crpix = _calculate_crpix(nside)
    wcs_.wcs.crpix = _crpix, _crpix
    # 2 element vector giving R.A. and DEC of reference pixel in degrees.
    wcs_.wcs.crval = ra0, dec0
    # 2 element string vector giving projection types.
    wcs_.wcs.ctype = 'RA---%s' % projection, 'DEC--%s' % projection
    # Miscellanea.
    wcs_.array_shape = nside, nside
    wcs_.wcs.equinox = 2000.
    wcs_.wcs.radesys = 'ICRS'
    wcs_.wcs.lonpole = 180.
    wcs_.wcs.latpole = 0.
    return wcs_


def wcs_to_sky_meshgrid(wcs_):
    """Convert a WCS object into a meshgrid of sky coordinates representing the
    (ra, dec) position of the pixel centers in the sky.

    This comes handy when evaluating a model over the sky grid defined by a
    WCS object.

    Arguments
    ---------
    wcs_ : astropy.wcs.WCS object
        The underlying WCS.
    """
    nx, ny = wcs_.array_shape
    # Define the positions of the bin centers in logical coordinates with
    # the Python convention (i.e., origin = 0).
    x = numpy.arange(0., nx, 1.)
    y = numpy.arange(0., ny, 1.)
    x, y = numpy.meshgrid(x, y, indexing='ij')
    return wcs_.wcs_pix2world(x, y, 0)


def wcs_digitize(wcs_, ra, dec, weights=None):
    """Digitize a set of sky coordinates over a given WCS.

    This returns an accumulated array with the same shape as the WCS, that can
    be readily plotted in the sky using the WCS itself.

    Arguments
    ---------
    wcs_ : astropy.wcs.WCS object
        The underlying WCS.

    ra : array_like
        The right ascension values.

    dec : array_like
        The declination values.

    weights : array_like, optional
        The weights to be passed to histogram2d.
    """
    # Hack for the scalar case.
    if isinstance(ra, numbers.Number):
        ra = numpy.array([ra], float)
    if isinstance(dec, numbers.Number):
        dec = numpy.array([dec], float)
    # At this point we must have ended up with numpy arrays of the same length.
    assert ra.shape == dec.shape
    nx, ny = wcs_.array_shape
    xbins = numpy.linspace(-0.5, nx - 0.5, nx + 1)
    ybins = numpy.linspace(-0.5, ny - 0.5, ny + 1)
    x, y = wcs_.wcs_world2pix(ra, dec, 0)
    data, _, _ = numpy.histogram2d(x, y, (xbins, ybins), weights=weights)
    # Mind to go along with the WCS the data need to be transposed.
    data = data.transpose()
    return data

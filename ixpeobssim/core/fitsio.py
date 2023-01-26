#!/usr/bin/env python
#
# Copyright (C) 2015--2020, the ixpeobssim team.
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

"""Facilities related to FITS I/O.
"""

from __future__ import print_function, division

import numpy
from matplotlib.patches import Circle
from matplotlib.ticker import LogFormatter, LogFormatterSciNotation, LogLocator

from astropy.io import fits
from astropy.wcs import wcs
from astropy.visualization.wcsaxes import WCSAxes
from astropy.visualization import simple_norm

from ixpeobssim.__version__ import TAG as IXPEOBSSIM_VERSION
from ixpeobssim.utils.time_ import current_datetime_string_utc
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.matplotlib_ import plt, plot_arrows
from ixpeobssim.utils.os_ import check_input_file
from ixpeobssim.utils.units_ import arcmin_to_degrees, degrees_to_arcsec
from ixpeobssim.utils.astro import square_sky_grid
from ixpeobssim.utils.matplotlib_ import draggable_colorbar



# pylint: disable=invalid-name, too-many-arguments, no-member, consider-using-f-string



# From https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node20.html
# Codes for the data type of binary table columns and/or for the
# data type of variables when reading or writing keywords or data:
#
#                               DATATYPE               TFORM CODE
#   #define TBIT          1  /*                            'X' */
#   #define TBYTE        11  /* 8-bit unsigned byte,       'B' */
#   #define TLOGICAL     14  /* logicals (int for keywords     */
#                            /*  and char for table cols   'L' */
#   #define TSTRING      16  /* ASCII string,              'A' */
#   #define TSHORT       21  /* signed short,              'I' */
#   #define TLONG        41  /* signed long,                   */
#   #define TLONGLONG    81  /* 64-bit long signed integer 'K' */
#   #define TFLOAT       42  /* single precision float,    'E' */
#   #define TDOUBLE      82  /* double precision float,    'D' */
#   #define TCOMPLEX     83  /* complex (pair of floats)   'C' */
#   #define TDBLCOMPLEX 163  /* double complex (2 doubles) 'M' */
#
#   The following data type codes are also supported by CFITSIO:
#   #define TINT         31  /* int                            */
#   #define TSBYTE       12  /* 8-bit signed byte,         'S' */
#   #define TUINT        30  /* unsigned int               'V' */
#   #define TUSHORT      20  /* unsigned short             'U'  */
#   #define TULONG       40  /* unsigned long                  */
#
#   The following data type code is only for use with fits_get_coltype
#   #define TINT32BIT    41  /* signed 32-bit int,         'J' */


FITS_TO_NUMPY_TYPE_DICT = {
    'E': numpy.float32,
    'D': numpy.float64,
    'I': numpy.int16,
    'J': numpy.int32
    }



def copy_hdu_list(hdu_list):
    """Create a copy in memory of a HDU list.

    Note that the HDUList class provides a copy() method, but this is a shallow
    copy that is still tied to the underlying file the HDU list is read from:
    https://github.com/astropy/astropy/issues/6395

    As a workaround we copy the HDUs one by one and assemble them back into a
    HDUList object.
    """
    return fits.HDUList([hdu.copy() for hdu in hdu_list])



def read_hdu_list_in_memory(file_path, extension='fits'):
    """Read a fits file and store the HDU list in memory.

    This is a convenience function that is meant to open a file within a context
    manager (so that it is properly closed after the fact) and create a copy
    of the content in memory. This is less than trivial, as astropy typically
    keeps some kind of pointer back to the original file when dealing with FITS
    files, and it is not obvious how to make sure that all the resources are
    properly closed at runtime.
    """
    check_input_file(file_path, extension)
    logger.info('Reading (in memory) %s...', file_path)
    with fits.open(file_path, lazy_load_hdus=False, mmap=False) as hdu_list:
        hdu_list = copy_hdu_list(hdu_list)
    return hdu_list



def find_column_index(hdu, col_name):
    """Return the index corresponding to a given column name within a binary table,
    and None when the column does not exists.

    Arguments
    ---------
    hdu : fits.BinTableHDU instance
        The underlying HDU

    col_name : str
        The name of the target column.
    """
    col_names = [col.name for col in hdu.data.columns]
    try:
        return col_names.index(col_name)
    except ValueError:
        logger.warning('Could not find column %s in extension %s.', col_name, hdu.name)
        return None



def set_tlbounds(hdu, col_name, min_, max_):
    """Set the proper TLMIN and TLMAX header keywords for a given column in a
    given HDU.
    """
    if not isinstance(hdu, fits.BinTableHDU):
        raise RuntimeError('HDU %s is not a binary table.' % hdu.name)
    index = find_column_index(hdu, col_name)
    if index is None:
        logger.error('Cannot set TLMIN/TLMAX for column %s in extension %s', col_name, hdu.name)
        return
    col_index = index + 1
    logger.info('Updating bounds for column %s in extension %s...', col_name, hdu.name)
    for key, value in zip(('TLMIN%d' % col_index, 'TLMAX%d' % col_index), (min_, max_)):
        logger.info('Setting %s -> %d', key, value)
        hdu.header.set(key, value)



class xHDUBase:

    """Base class for FITS HDU.
    """

    def add_comment(self, comment):
        """Add a comment to the table header.
        """
        self.header['COMMENT'] = comment

    def add_keyword(self, key, value, comment=''):
        """Add a keyword to the table header.
        """
        self.header.set(key, value, comment)

    def set_keyword(self, key, value):
        """Set a keyword into the table header.
        """
        self.header.set(key, value)

    def set_keyword_comment(self, keyword, comment):
        """Set the comment for a header keyword.
        """
        self.header.comments[keyword] = comment

    def setup_header(self, keywords=None, comments=None):
        """Update the table header with arbitrary additional information.
        """
        if keywords is None:
            keywords = []
        if comments is None:
            comments = []
        for item in keywords:
            if len(item) == 3:
                key, value, comment = item
            elif len(item) == 2:
                key, value = item
                comment = ''
            self.add_keyword(key, value, comment)
        for comment in comments:
            self.add_comment(comment)

    def __str__(self):
        """String formatting.
        """
        return repr(self.header)



class xPrimaryHDU(fits.PrimaryHDU, xHDUBase):

    """Class describing a primary HDU to be written in a FITS file.

    This is initializing a standard astropy.io.fits.PrimaryHDU object and
    adding the creator and creation time fields.

    Arguments
    ---------
    creator : str
       The application that created the header (defaults to `ixpeobssim`, and
       the ixpeobssim tag is automatically added.)
    """

    HEADER_KEYWORDS = []
    HEADER_COMMENTS = []

    def __init__(self, data=None, header=None, creator='ixpeobssim', keywords=None,
                 comments=None):
        """Constructor.
        """
        if keywords is None:
            keywords = []
        if comments is None:
            comments = []
        fits.PrimaryHDU.__init__(self, data, header)
        self.header.set('CREATOR', '%s ver. %s' % (creator, IXPEOBSSIM_VERSION),
                        's/w task which wrote this dataset')
        self.header.set('DATE', current_datetime_string_utc(), 'file creation date')
        self.setup_header(self.HEADER_KEYWORDS + keywords,
                          self.HEADER_COMMENTS + comments)



class xBinTableHDUBase(fits.BinTableHDU, xHDUBase):

    """Binary table HDU class.

    This is a small wrapper around a standard binary table to facilitate
    customizations.
    """

    NAME = None
    HEADER_KEYWORDS = []
    HEADER_COMMENTS = []
    DATA_SPECS = []
    DATA_KWARGS = None

    def __init__(self, data=None, keywords=None, comments=None):
        """
        """
        if keywords is None:
            keywords = []
        if comments is None:
            comments = []
        if data is not None:
            # The data can be passed in the form of a dictionary, in which the
            # columns are automatically extracted.
            if isinstance(data, dict):
                data = [data[name] for name in self.spec_names()]
            if len(data) != len(self.DATA_SPECS):
                logger.error(data)
                logger.error(self.DATA_SPECS)
                abort('Data-specs length mismatch (%d vs. %s)' % (len(data), len(self.DATA_SPECS)))
        cols = []
        _kwcomments = {}
        for i, item in enumerate(self.DATA_SPECS):
            if len(item) == 4:
                name, format_, units, comment = item
                _kwcomments[name] = comment
            elif len(item) == 3:
                name, format_, units = item
            elif len(item) == 2:
                name, format_ = item
                units = None
            # New switch to support the necessary keyword arguments for the
            # X and Y columns
            if self.DATA_KWARGS is not None:
                _kwargs = self.DATA_KWARGS.get(name, {})
            else:
                _kwargs = {}
            if len(_kwargs):
                logger.info('Passing %s keyword argument(s) to column %s...', _kwargs, name)
            if data is not None:
                col = fits.Column(name, format_, units, array=data[i], **_kwargs)
            else:
                col = fits.Column(name, format_, units, **_kwargs)
            cols.append(col)
        data = fits.FITS_rec.from_columns(cols)
        fits.BinTableHDU.__init__(self, data)
        # Set the extension name, if necessary.
        if self.NAME is not None:
            self.set_ext_name(self.NAME)
        # Add the additional keywords and comments to the header.
        self.setup_header(self.HEADER_KEYWORDS + keywords,
                          self.HEADER_COMMENTS + comments)
        # And we need to loop one more time to take care of the comments on the
        # columns, if any (could not find a way to do this on the columns
        # directly).
        for i, col in enumerate(self.columns):
            if col.name in _kwcomments:
                comment = _kwcomments[col.name]
                self.set_keyword_comment('TTYPE%d' % (i + 1), comment)

    @classmethod
    def spec_names(cls):
        """Return the name of the data fields specified in the SPEC class
        member.
        """
        return [item[0] for item in cls.DATA_SPECS]

    @classmethod
    def spec_names_and_types(cls):
        """Return the name of the data fields specified in the SPEC class
        member.
        """
        return [item[0:2] for item in cls.DATA_SPECS]

    def set_ext_name(self, name):
        """Set the extension name for the binary table.
        """
        self.add_keyword('EXTNAME', name, 'Name of this binary table extension')

    def __str__(self):
        """String formatting.
        """
        return '%s\n%s' % (repr(self.header), self.data)



class xFITSImageBase:

    """Base class describing a FITS image.

    Arguments
    ---------
    file_path : string
        The path to the FITS file containing the image.
    """

    DEFAULT_RECT = (0.1, 0.1, 0.8, 0.8)

    def __init__(self, file_path):
        """Constructor.
        """
        self.hdu_list = read_hdu_list_in_memory(file_path)
        self.hdu_list.info()
        self.primary_hdu = self.hdu_list['PRIMARY'].copy()
        self.wcs = wcs.WCS(self.primary_hdu.header)
        self.data = self.primary_hdu.data

    def __iadd__(self, other):
        """Support for image addition.

        This comes handy, e.g., when loading binned maps in FITS format for the
        three Detector units and displaying the combined map.
        """
        assert isinstance(other, self.__class__)
        assert self.data.shape == other.data.shape
        self.data += other.data
        return self

    def get(self, keyword, default=None):
        """Retrieve the value of a primary header keyword.
        """
        return self.primary_hdu.header.get(keyword, default)

    def shape(self):
        """Return the shape of the image.
        """
        return self.data.shape

    def crval1(self):
        """Return the value of the CRVAL1 header keyword.
        """
        return self.get('CRVAL1')

    def crval2(self):
        """Return the value of the CRVAL2 header keyword.
        """
        return self.get('CRVAL2')

    def crpix1(self):
        """Return the value of the CRPIX1 header keyword.
        """
        return self.get('CRPIX1')

    def crpix2(self):
        """Return the value of the CRPIX2 header keyword.
        """
        return self.get('CRPIX2')

    def cdelt1(self):
        """Return the value of the CRVAL1 header keyword.
        """
        return self.get('CDELT1')

    def cdelt2(self):
        """Return the value of the CRVAL2 header keyword.
        """
        return self.get('CDELT2')

    def delta(self):
        """Return the pixel increment in the ra and dec coorinates.
        """
        return numpy.array([self.cdelt1(), self.cdelt2()])

    def inner_radius(self):
        """Return the inner radius of the image in degrees.
        """
        return 0.5 * abs(self.delta() * self.shape()).min()

    def sky_bounding_box(self):
        """Return the bounding box of the FITS image in sky coordinates, i.e.,
        a 4-element tuple of the form (ra_min, dec_min, ra_max, dec_max).
        """
        nx, ny = self.shape()
        ra_min, dec_min = self.wcs.wcs_pix2world(-0.5, -0.5, 0)
        ra_max, dec_max = self.wcs.wcs_pix2world(nx - 0.5, ny - 0.5, 0)
        return ra_min, dec_min, ra_max, dec_max

    def center(self):
        """Return the sky coordinates of the image center.
        """
        nx, ny = self.shape()
        x, y = 0.5 * nx - 1., 0.5 * ny - 1.
        return self.wcs.wcs_pix2world(x, y, 0)

    def __call__(self, i, j):
        """Return the value of the underlying map for a given pixel.

        Note the addressing is in FITS space, i.e., one pixel unit is subtracted
        in both dimensions to address the underlying numpy array.

        Warning
        -------
        I am not sure this is even used, and given the ambiguity we have with
        the wcs origin in general, I would be in favor of deprecating this.
        """
        return self.data[i - 1][j - 1]

    @staticmethod
    def make_plot(data, wcs_, slices=None, zlabel=None, stretch='linear',
                  ticks=None, colorbar_format=None, **kwargs):
        """Underlying plotting routine.

        This is implemented as a staticmethod in such a way it can be called
        from external methods dealing with multi-dimensional data arrays, which
        are not readily supported in the base class. This is achieved by
        passing the slices argument, aling with the proper data slice.

        Arguments
        ---------
        data : array_like
            The underlying data to be plotted.

        wcs_ : astropy.wcs object
            The WCS object defining the transformation from pixel to sky coordinates.

        slices : tuple, optional
            Verbatim from the astropy.visualization documentation:
            For WCS transformations with more than two dimensions, we need to
            choose which dimensions are being shown in the 2D image. The slice
            should contain one x entry, one y entry, and the rest of the values
            should be integers indicating the slice through the data. The order
            of the items in the slice should be the same as the order of the
            dimensions in the WCS, and the opposite of the order of the
            dimensions in Numpy. For example, (50, 'x', 'y') means that the
            first WCS dimension (last Numpy dimension) will be sliced at an
            index of 50, the second WCS and Numpy dimension will be shown on
            the x axis, and the final WCS dimension (first Numpy dimension)
            will be shown on the y-axis (and therefore the data will be plotted
            using data[:, :, 50].transpose())

        zlabel : str
            The label for the z axis (e.g., the colorbar).

        stretch : {‘linear’, ‘sqrt’, ‘power’, log’, ‘asinh’}, optional
            The stretch function to apply to the image.

        ticks : array_like or ticker, optional
            The explicit setting for the colorbar ticks.

        colorbar_format : {None, 'log', 'scilog'}, optional
            Optional string for the colorbar formnatting.

        kwargs : dict
            All the keyword arguments passed to plt.imshow().
        """
        # Set sensible overall default.
        kwargs.setdefault('cmap', 'afmhot')
        kwargs.setdefault('aspect', 'equal')
        # Handle vmin and vmax. This needs special care, as since matplotlib
        # 3.3 vmin and vmax shound not be passed as kwargs to functions using
        # colormappings when norm is used, as in our case.
        vmin = kwargs.pop('vmin', data.min())
        vmax = kwargs.pop('vmax', data.max())
        # If the stretch is logarithmic we do some magic to make sure we get
        # sensible default settings.
        if stretch == 'log':
            # For a log stretch we have to make sure that vmin > 0.
            if vmin == 0.:
                try:
                    # Note this might fail if the image has no positive values
                    # (e.g., if it is identically zero.)
                    vmin = data[data > 0.].min()
                except ValueError:
                    # In which case we literally have to make something up.
                    vmin = 0.1
                    if vmax <= vmin:
                        vmax = 1.
            # If we're not setting the formatter for the color bar, use a
            # logarithmic formatted in scientific notation.
            if colorbar_format is None:
                colorbar_format = 'scilog'
            # If we're not setting the locator for the color bar, use a
            # plain logarithmic locator.
            if ticks is None:
                ticks = LogLocator()
        # Set up the image normalization. Note we adjust the log_a argument
        # simple_norm() in such a way the powers of ten are equispaced on the screen.
        log_a = vmax / vmin
        norm = simple_norm(data, stretch, clip=True, log_a=log_a, min_cut=vmin, max_cut=vmax)
        # Set up the formatter for the color bar.
        if colorbar_format == 'scilog':
            colorbar_format = LogFormatterSciNotation(labelOnlyBase=True)
        elif colorbar_format == 'log':
            colorbar_format = LogFormatter(labelOnlyBase=True)
        # We're ready to move on with the actual figure.
        fig = plt.gcf()
        ax = WCSAxes(fig, xFITSImageBase.DEFAULT_RECT, wcs=wcs_, slices=slices)
        fig.add_axes(ax)
        im = ax.imshow(data, norm=norm, **kwargs)
        ax.set_xlabel('Right Ascension (J2000)')
        ax.set_ylabel('Declination (J2000)')
        ax.grid(color='gray')
        colorbar = draggable_colorbar(im, format=colorbar_format, ticks=ticks, pad=0.05)
        if zlabel is not None:
            colorbar.set_label(zlabel)
        return fig

    def plot(self, zlabel=None, stretch='linear', ticks=None,
             colorbar_format=None, **kwargs):
        """Plot the underlying FITS image.

        Since ixpeobssim version 8.9.0 this not using aplpy under the hood
        anymore, see https://bitbucket.org/ixpesw/ixpeobssim/issues/272

        Arguments
        ---------
        zlabel : string
            The text label for the colorbar (use None for no colorbar)

        **kwargs
            All the keyword arguments passed to plt.imshow()
        """
        return self.make_plot(self.data, self.wcs, None, zlabel, stretch, ticks,
                              colorbar_format, **kwargs)

    @staticmethod
    def add_label(text, x=0.1, y=0.9, **kwargs):
        """Add a label to an image.

        This is a shortcut to have all the formatting defined.

        Arguments
        ---------
        text : string
            The label text

        x : float
            The x position of the label in relative coordinates

        y : float
            The y position of the label in relative coordinates

        **kwargs
            All the keyword arguments passed to plt.text()
        """
        kwargs.setdefault('color', 'white')
        kwargs.setdefault('size', 'large')
        kwargs.setdefault('ha', 'left')
        kwargs.setdefault('transform', plt.gca().transAxes)
        plt.text(x, y, text, **kwargs)

    def recenter(self, ra, dec, radius):
        """Recenter the image.

        Arguments
        ---------
        ra : float
            The RA position of the circle center in decimal degrees.

        dec : float
            The DEC position of the circle center in decimal degrees.

        radius : float
            The radius of the region in arcseconds.
        """
        x, y = self.wcs.wcs_world2pix(ra, dec, 0)
        dx = radius / degrees_to_arcsec(abs(self.cdelt1()))
        dy = radius / degrees_to_arcsec(abs(self.cdelt2()))
        plt.gca().axis([x - dx, x + dx, y - dy, y + dy])

    def add_circle(self, x, y, radius, mode='radec', **kwargs):
        """Add a circle to the image.

        Note that the x and y coordinates can be either RA and Dec in decimal
        degrees (if mode='radec') or relative to canvas coordinates
        (if mode='canvas').

        Arguments
        ---------
        x : float
            The x position of the circle center (Right Ascension in decimal
            degrees if relative=False, relative position in canvas coordinates
            if relative=True).

        y : float
            The y position of the circle center (Declination in decimal degrees
            if relative=False, relative position in canvas coordinates if
            relative=True).

        radius : float
            The radius of the circle in arcseconds.

        mode : {'radec', 'canvas'}
            Optional flag defining how the center of the circle is expresses:
            in right ascension and declination ('radec') or in relative
            canvas coordinates ('canvas').

        kwargs : dict
            All the keyword arguments passed to the matplotlib Circle patch.
        """
        kwargs.setdefault('color', 'white')
        kwargs.setdefault('fill', None)
        assert mode in ('radec', 'canvas')
        if mode == 'radec':
            x, y = self.wcs.wcs_world2pix(x, y, 0)
            radius /= degrees_to_arcsec(abs(self.cdelt2()))
        elif mode == 'canvas':
            kwargs.setdefault('transform', plt.gca().transAxes)
            radius /= degrees_to_arcsec(self.inner_radius())
        circle = Circle((x, y), radius, **kwargs)
        plt.gca().add_patch(circle)

    def add_circle_radec(self, ra, dec, radius, **kwargs):
        """Convenience function, see add_circle().
        """
        self.add_circle(ra, dec, radius, mode='radec', **kwargs)

    def add_circle_canvas(self, ra, dec, radius, **kwargs):
        """Convenience function, see add_circle().
        """
        self.add_circle(ra, dec, radius, mode='canvas', **kwargs)

    def square_sky_grid(self, nside, ra0=None, dec0=None, half_size=None):
        """Create a square, regular grid in ra and dec.

        This is essentially an overloaded version of the
        utils.astro.square_sky_grid() method, where the arguments are inferred,
        when possible, from the underlying image structure.
        """
        if ra0 is None:
            ra0 = self.crval1()
        if dec0 is None:
            dec0 = self.crval2()
        if half_size is None:
            half_size = self.inner_radius()
        else:
            half_size = arcmin_to_degrees(half_size)
        return square_sky_grid(nside, (ra0, dec0), half_size)

    def plot_arrows(self, nside, model, threshold=0., **kwargs):
        """Overlay arrows on an image.
        """
        plot_arrows(self.square_sky_grid(nside), model, threshold, **kwargs)

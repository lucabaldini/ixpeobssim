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


"""Event display facilities.

This module provides a simple, top-level interface to track images in IXPE
level-1 files.
"""

from dataclasses import dataclass

from astropy.io import fits
import numpy
import matplotlib

from matplotlib.patches import RegularPolygon
from matplotlib.collections import PatchCollection

from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt


XPOL_LAYOUT = 'ODD_R'
XPOL_SIZE = (300, 352)
XPOL_VERTICAL_PADDING = 10
XPOL_HORIZONTAL_PADDING = 8
XPOL_PITCH = 0.05 # mm

# pylint: disable = invalid-name



@dataclass
class xRegionOfInterest:

    """Class describing a region of interest (ROI).

    A region of interest is the datum of the logical coorinates of its two
    extreme corners, in the order (min_col, max_col, min_row, max_row).
    """

    min_col : int
    max_col : int
    min_row : int
    max_row : int

    def __post_init__(self):
        """Overloaded dataclass method.
        """
        self.num_cols = self.max_col - self.min_col + 1
        self.num_rows = self.max_row - self.min_row + 1
        self.size =  self.num_cols * self.num_rows
        self.shape = (self.num_rows, self.num_cols)

    def at_border(self):
        """Return True if the ROI is on the border for a given chip_size.
        """
        num_cols, num_rows = XPOL_SIZE
        return self.min_col == 0 or self.max_col == num_cols - 1 or\
            self.min_row == 0 or self.max_row == num_rows - 1

    def column_indices(self):
        """Return an array with all the valid column indices.
        """
        return numpy.arange(self.min_col, self.max_col + 1)

    def row_indices(self):
        """Return an array with all the valid row indices.
        """
        return numpy.arange(self.min_row, self.max_row + 1)

    def serial_readout_coordinates(self):
        """Return two one-dimensional arrays containing the column and row
        indices, respectively, in order of serial readout of the ROI.

        Example
        -------
        >>> col, row = xRegionOfInterest(0, 4, 0, 3).serial_readout_coordinates()
        >>> print(col)
        >>> [0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
        >>> print(row)
        >>> [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]
        """
        col = numpy.tile(self.column_indices(), self.num_rows)
        row = numpy.repeat(self.row_indices(), self.num_cols)
        return col, row

    def coordinates_in_rot(self, col, row):
        """Return a boolean mask indicaing whether elements of the (col, row)
        arrays lie into the ROT area.
        """
        return numpy.logical_and.reduce((
            col >= self.min_col + XPOL_HORIZONTAL_PADDING,
            col <= self.max_col - XPOL_HORIZONTAL_PADDING,
            row >= self.min_row + XPOL_VERTICAL_PADDING,
            row <= self.max_row - XPOL_VERTICAL_PADDING
        ))



@dataclass
class Recon:

    """Container class encapsulating the event reconstruction.
    """

    absorption_point : tuple[float, float]
    barycenter : tuple[float, float]
    track_direction : float
    length : float
    width : float

    @staticmethod
    def _annotate_point(x, y, text, xoffset=125, yoffset=75, color='black'):
        """Small convenience function to annotate a point.
        """
        plt.plot(x, y, 'o', color=color)
        arrowprops=dict(arrowstyle='-', connectionstyle='angle3', color=color)
        kwargs = dict(xycoords='data', textcoords='offset points', arrowprops=arrowprops,
            backgroundcolor='white', color=color, ha='center')
        plt.gca().annotate(text, xy=(x, y), xytext=(xoffset, yoffset), **kwargs)

    def draw_absorption_point(self):
        """Draw the reconstructed absorption point and track direction.
        """
        xoffset = 125 if self.absorption_point[0] >= self.barycenter[0] else -125
        self._annotate_point(*self.absorption_point, 'Absorption point', xoffset)

    def draw_barycenter(self):
        """Draw the reconstructed absorption point and track direction.
        """
        xoffset = 125 if self.absorption_point[0] < self.barycenter[0] else -125
        self._annotate_point(*self.barycenter, 'Barycenter', xoffset, color='#777')

    def draw_track_direction(self, line_width=1., length_ratio=0.5):
        """Draw the track direction.
        """
        x0, y0 = self.absorption_point
        dist = max(0.5, self.length)
        dx = dist * numpy.cos(self.track_direction)
        dy = dist * numpy.sin(self.track_direction)
        if dist > 1:
            length_ratio /= dist
        plt.plot([x0, x0 - length_ratio * dx], [y0, y0 - length_ratio * dy],
            lw=line_width, color='black', ls='dashed')
        plt.gca().annotate('', xy=(x0 + dx, y0 + dy), xytext=(x0, y0),
            arrowprops=dict(arrowstyle='->, head_width=0.25, head_length=0.5', lw=line_width))


@dataclass
class xL1Event(xRegionOfInterest):

    """A fully fledged event.

    This is building up on the core logic encapsulated in the xRegionOfInterest
    class, and is adding all the necessary event information on top of that, i.e.,
    the pha values, and the time and error information.

    Note the dataclass field are largely mapped over the columns of the
    corresponding EVENTS extension in the underlying FITS files.
    """

    pha : numpy.array
    trigger_id : int = 0
    seconds : int = 0
    microseconds : int = 0
    timestamp : float = 0.
    livetime : int = 0
    error_summary : int = 0
    du_status : int = 0
    recon : Recon = None

    def __post_init__(self):
        """Post-init hook implementation.
        """
        super().__post_init__()
        self.pha = self.pha.reshape(self.shape)

    def highest_pixel(self, absolute=True):
        """Return the coordinates (col, row) of the highest pixel.

        Arguments
        ---------
        absolute : bool
            If true, the absolute coordinates (i.e., those referring to the readout
            chip) are returned; otherwise the coordinates are intended relative
            to the readout window (i.e., they can be used to index the pha array).
        """
        # Note col and row are swapped, here, due to how the numpy array are indexed.
        # pylint: disable = unbalanced-tuple-unpacking
        row, col = numpy.unravel_index(numpy.argmax(self.pha), self.pha.shape)
        if absolute:
            col += self.min_col
            row += self.min_row
        return col, row



# pylint: disable = too-few-public-methods
class xL1EventFile:

    """Simple interface to a level-1 file in FITS format.

    Arguments
    ---------
    file_path : str
        The path to the input event file.

    padding : Padding instance
        The ROI padding for the event file.
    """

    EVT_EXT_NAME = 'EVENTS'
    EVT_COL_NAMES = ('MIN_CHIPX', 'MAX_CHIPX', 'MIN_CHIPY', 'MAX_CHIPY', 'PIX_PHAS',
        'TRG_ID', 'SEC', 'MICROSEC', 'TIME', 'LIVETIME', 'ERR_SUM', 'DU_STATUS')

    def __init__(self, file_path):
        """Constructor.
        """
        # pylint: disable = no-member
        logger.info('Opening input file %s...', file_path)
        self.hdu_list = fits.open(file_path)
        self.__index = -1
        self.__num_events = len(self.hdu_list[self.EVT_EXT_NAME].data)
        # Cache the event MET values, since we typically need to bisect downstream.
        self.__met_values = self.hdu_list[self.EVT_EXT_NAME].data['TIME']
        logger.info('Done, %d event(s) found.', self.__num_events)
        logger.warning('Mind that indexing the events might take some time...')

    def zero_sup_threshold(self):
        """Return the zero-suppression threshold, determined from the file header.
        """
        return self.hdu_list[self.EVT_EXT_NAME].header['ZSUPTHR']

    def value(self, col_name):
        """Return the value of a given column for a given extension for the current event.
        """
        return self.hdu_list[self.EVT_EXT_NAME].data[col_name][self.__index]

    def _recon(self):
        """Retrieve the reconstructed quantities.
        """
        absorption_point = self.value('ABSX'), self.value('ABSY')
        barycenter = self.value('BARX'), self.value('BARY')
        track_direction = self.value('DETPHI2')
        length = 4. * numpy.sqrt(self.value('TRK_M2L'))
        width = 4. * numpy.sqrt(self.value('TRK_M2T'))
        return Recon(absorption_point, barycenter, track_direction, length, width)

    def __getitem__(self, event_number):
        """Overloaded slicing hook.

        This returns a fully fledged Event object for a given event number.
        """
        self.__index = event_number
        args = [self.value(col_name) for col_name in self.EVT_COL_NAMES]
        return xL1Event(*args, self._recon())

    def bisect_met(self, met):
        """Retrieve a specific event by its mission elapsed time.

        Internally this is using a binary search on the time column, and in
        general it can be assumed that this O(log(N)) in complexity.
        """
        event_number = numpy.searchsorted(self.__met_values, met)
        event = self[event_number]
        delta_time = event.timestamp - met
        if abs(delta_time) > 1.e-6:
            logger.warning('Bisected MET is %.6f s (target: %.6f s, difference %.6f s)',
                event.timestamp, met, delta_time)
        self.__index = event_number
        return event

    def __iter__(self):
        """Iterator protocol implementation.
        """
        return self

    def __next__(self):
        """Iterator protocol implementation.
        """
        self.__index += 1
        if self.__index >= self.__num_events:
            raise StopIteration
        return self[self.__index]



class xHexagonCollection(PatchCollection):

    """Collection of native hexagin patches.

    Arguments
    ---------
    x : array_like
        The x coordinates of the hexagon centers.

    y : array_like
        The y coordinates of the hexagon centers.

    radius : float
        The hexagon apothem.

    orientation: float
        The hexagon orientation in radians---zero means pointy topped.

    kwargs
        The keyword arguments to be passed to the PatchCollection constructor.
    """

    DEFAULT_EDGE_COLOR = '#CCC'

    def __init__(self, x, y, radius=XPOL_PITCH, orientation=0., **kwargs):
        """Constructor.
        """
        # pylint: disable = invalid-name
        self.x = x
        self.y = y
        kwargs.setdefault('edgecolor', self.DEFAULT_EDGE_COLOR)
        kwargs.setdefault('facecolor', 'none')
        kwargs.setdefault('linewidth', 1.2)
        patches = [RegularPolygon(xy, 6, radius, orientation) for xy in zip(x, y)]
        PatchCollection.__init__(self, patches, match_original=False, **kwargs)



class xHexagonalGrid:

    """Generic hexagonal grid.

    Arguments
    ---------
    num_cols : int
        The number of columns in the grid

    num_rows : int
        The number of rows in the grid

    pitch : float
        The grid pitch in mm.
    """

    # pylint: disable = too-many-instance-attributes
    def __init__(self, num_cols, num_rows, pitch=XPOL_PITCH, **kwargs):
        """Constructor.
        """
        self.num_cols = num_cols
        self.num_rows = num_rows
        self.pitch = pitch
        self.color_map = matplotlib.cm.get_cmap(kwargs.get('cmap_name', 'Reds')).copy()
        self.color_map_offset = kwargs.get('cmap_offset', 0)
        self.color_map.set_under('white')
        self.secondary_pitch = 0.5 * numpy.sqrt(3.) * self.pitch
        self.xoffset = 0.5 * (self.num_cols - 1.5) * self.pitch
        self.yoffset = 0.5 * (self.num_rows - 1) * self.secondary_pitch

    def pixel_to_world(self, col, row):
        """Transform pixel to world coordinates.

        Arguments
        ---------
        col : array_like
            The input column number(s).

        row : array_like
            The input row number(s).
        """
        # pylint: disable = invalid-name
        x = (col -0.5 * (row & 1)) * self.pitch - self.xoffset
        y = self.yoffset - row * self.secondary_pitch
        return x, y

    def pha_to_colors(self, pha, zero_sup_threshold=None):
        """Convert the pha values to colors for display purposes.
        """
        values = pha.flatten()
        values += self.color_map_offset
        if zero_sup_threshold is not None:
            values[values <= zero_sup_threshold + self.color_map_offset] = -1.
        values = values / float(values.max())
        return self.color_map(values)

    # pylint: disable = too-many-arguments, too-many-locals
    def draw_roi(self, roi, offset=(0., 0.), indices=True, padding=True, **kwargs):
        """Draw a specific ROI of the parent grid.
        """
        # pylint: disable = invalid-name
        # Calculate the coordinates of the pixel centers and build the basic
        # hexagon collection.
        col, row = roi.serial_readout_coordinates()
        dx, dy = offset
        x, y = self.pixel_to_world(col, row)
        collection = xHexagonCollection(x + dx, y + dy, 0.5 * self.pitch, **kwargs)
        # If the padding is defined, we want to distinguish the different regions
        # by the pixel edge color.
        if padding:
            color = numpy.full(col.shape, '#555')
            color[~roi.coordinates_in_rot(col, row)] = xHexagonCollection.DEFAULT_EDGE_COLOR
            collection.set_edgecolor(color)
        plt.gca().add_collection(collection)
        # And if we want the indices, we add appropriate text patches.
        if indices:
            font_size = 'x-small'
            cols, rows = roi.column_indices(), roi.row_indices()
            first_row = numpy.full(cols.shape, roi.min_row)
            first_col = numpy.full(rows.shape, roi.min_col)
            fmt = dict(fontsize=font_size, ha='center', va='bottom', rotation=60.)
            for x, y, col in zip(*self.pixel_to_world(cols, first_row), cols):
                plt.text(x + dx, y + dy + self.secondary_pitch, f'{col}', **fmt)
            fmt = dict(fontsize=font_size, ha='right', va='center', rotation=0.)
            for x, y, row in zip(*self.pixel_to_world(first_col, rows), rows):
                plt.text(x + dx - self.pitch, y + dy, f'{row}', **fmt)
        return collection

    @staticmethod
    def brightness(color):
        """Quick and dirty proxy for the brighness of a given array of colors.

        See https://stackoverflow.com/questions/9733288
        and also
        https://stackoverflow.com/questions/30820962
        for how to split in columns the array of colors.
        """
        # pylint: disable = invalid-name
        r, g, b, _ = color.T
        return (299 * r + 587 * g + 114 * b) / 1000

    def draw_event(self, event, offset=(0., 0.), canvas_side=2.0, indices=True,
                   padding=True, zero_sup_threshold=None, values=False, **kwargs):
        """Draw an actual event int the parent hexagonal grid.

        This is taking over where the draw_roi() hook left, and adding the
        event part.
        """
        # pylint: disable = invalid-name
        collection = self.draw_roi(event, offset, indices, padding, **kwargs)
        face_color = self.pha_to_colors(event.pha, zero_sup_threshold)
        collection.set_facecolor(face_color)
        if values:
            # Draw the pixel values---note that we use black or white for the text
            # color depending on the brightness of the pixel.
            black = numpy.array([0., 0., 0., 1.])
            white = numpy.array([1., 1., 1., 1.])
            text_color = numpy.tile(black, len(face_color)).reshape(face_color.shape)
            text_color[self.brightness(face_color) < 0.5] = white
            fmt = dict(ha='center', va='center', fontsize='xx-small')
            for x, y, value, color in zip(collection.x, collection.y, event.pha.flatten(),
                text_color):
                if value > zero_sup_threshold:
                    plt.text(x, y, f'{value}', color=color, **fmt)
        if canvas_side is not None:
            x0, y0 = event.recon.barycenter
            pad = 0.5 * canvas_side
            plt.gca().set_xlim(x0 - pad, x0 + pad)
            plt.gca().set_ylim(y0 - pad, y0 + pad)
        return collection

    @staticmethod
    def show_display():
        """Convenience function to setup the matplotlib canvas for an event display.
        """
        plt.gca().set_aspect('equal')
        plt.axis('off')
        logger.info('Showing event display, close the window to move to the next one...')
        plt.show()



class xXpolGrid(xHexagonalGrid):

    """XPOL grid.
    """

    def __init__(self, **kwargs):
        """Constructor.
        """
        super().__init__(*XPOL_SIZE, XPOL_PITCH, **kwargs)

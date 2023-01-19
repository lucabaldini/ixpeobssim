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

from __future__ import annotations

from dataclasses import dataclass
import os

from astropy.io import fits
import numpy
import matplotlib

from matplotlib.patches import RegularPolygon
from matplotlib.collections import PatchCollection

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.evt.clustering import region_query_factory
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.matplotlib_ import plt, xStatBox, xTextCard


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

    def serial_readout_indices(self):
        """Return a two-dimensional array containing the readout index for
        each pixel of the ROI.

        The ASIC serial readout starts from the top-left corner and proceeds
        one row at a time, i.e., for a toy 5 x 5 window starting at <0, 0> the
        readout indices look like
        +--------------------
        |   0   1   2   3   4
        |
        |   5   6   7   8   9
        |
        |  10  11  12  13  14
        |
        |  15  16  17  18  19
        |
        |  20  21  22  23  24

        It is worth noting that, provided that one loops over the row indices
        first and column indices last, the readout index can be determined by
        just accumulating a counter. This function is useful as it provides the
        right answer for any position in the event window, no matter what the
        loop order is."""
        return numpy.arange(self.size).reshape(self.shape)

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

    def coordinates_in_roi(self, col, row):
        """Return a boolean mask indicaing whether elements of the (col, row)
        arrays lie into the ROT area.
        """
        return numpy.logical_and.reduce((
            col >= self.min_col,
            col <= self.max_col,
            row >= self.min_row,
            row <= self.max_row
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
        plt.plot(x, y, 'o', color=color, markersize=9.)
        arrowprops=dict(arrowstyle='-', connectionstyle='angle3', color=color)
        kwargs = dict(xycoords='data', textcoords='offset points', arrowprops=arrowprops,
            backgroundcolor='white', color=color, ha='center',
            bbox=dict(boxstyle='square,pad=0.', fc='white', ec='none'))
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
            arrowprops=dict(arrowstyle='->, head_width=0.4, head_length=0.65',
            lw=line_width))


@dataclass
class xL1Event(xRegionOfInterest):

    """A fully fledged event.

    This is building up on the core logic encapsulated in the xRegionOfInterest
    class, and is adding all the necessary event information on top of that, i.e.,
    the pha values, and the time and error information.

    Note the dataclass field are largely mapped over the columns of the
    corresponding EVENTS extension in the underlying FITS files.
    """

    _DIAGNOSTIC_DU_STATUS_BIT = 1
    _DIAGNOSTIC_OFFSET = 256

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
        # Handle diagnostic events.
        if (self.du_status >> self._DIAGNOSTIC_DU_STATUS_BIT) & 0x1:
            self.pha -= self._DIAGNOSTIC_OFFSET
        self.pha = self.pha.reshape(self.shape)
        self.cluster_id = numpy.zeros(self.shape, dtype=int)

    def run_clustering(self, engine):
        """Run the clustering on the track image.
        """
        region_query = region_query_factory(self)
        cluster_id = self.cluster_id.flatten()
        engine.run(self.pha.flatten(), cluster_id, region_query)
        self.cluster_id = cluster_id.reshape(self.shape)

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
        try:
            return self.hdu_list[self.EVT_EXT_NAME].data[col_name][self.__index]
        except IndexError as e:
            logger.warning(e)
            return None

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

    def roi_center(self, roi):
        """Return the world coordinates of the physical center of a given ROI.

        Note the small offsets that we apply serves the purpose of taking into
        account the fact that pixels are staggered in one direction.
        """
        col = (roi.min_col + roi.max_col + 1) // 2
        row = (roi.min_row + roi.max_row + 1) // 2
        x, y = self.pixel_to_world(col, row)
        x -= 0.75 * self.pitch
        y += 0.5 * self.secondary_pitch
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

    def default_roi_side(self, roi, min_side, pad=0.1):
        """Return the default physical size of the canvas necessary to fully
        contain a given ROI.
        """
        return pad + max(self.pitch * roi.num_cols, self.secondary_pitch * roi.num_rows,
            min_side)

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

    def draw_event(self, event, num_clusters=1, offset=(0., 0.), min_canvas_side=2.0, indices=True,
                   padding=True, zero_sup_threshold=None, values=False, **kwargs):
        """Draw an actual event int the parent hexagonal grid.

        This is taking over where the draw_roi() hook left, and adding the
        event part.
        """
        # pylint: disable = invalid-name
        # Create a copy of the PHA vector and set to zero all the pixels not
        # belonging to the target cluster.
        pha = event.pha.copy()
        # Note that negative numbers are used for orphan pixels, and if we want to
        # plot the first n clusters, we want the cluster id associated to the
        # pixel to be between 0 and n - 1.
        mask = numpy.logical_and(event.cluster_id >= 0, event.cluster_id < num_clusters)
        pha[numpy.logical_not(mask)] = 0
        # We're good to go!
        collection = self.draw_roi(event, offset, indices, padding, **kwargs)
        face_color = self.pha_to_colors(pha, zero_sup_threshold)
        collection.set_facecolor(face_color)
        if values:
            # Draw the pixel values---note that we use black or white for the text
            # color depending on the brightness of the pixel.
            black = numpy.array([0., 0., 0., 1.])
            white = numpy.array([1., 1., 1., 1.])
            text_color = numpy.tile(black, len(face_color)).reshape(face_color.shape)
            text_color[self.brightness(face_color) < 0.5] = white
            fmt = dict(ha='center', va='center', fontsize='xx-small')
            for x, y, value, color in zip(collection.x, collection.y, pha.flatten(),
                text_color):
                if value > zero_sup_threshold:
                    plt.text(x, y, f'{value}', color=color, **fmt)
        canvas_side = self.default_roi_side(event, min_canvas_side)
        # We want to center the display of the geometrical center of the ROI.
        x0, y0 = self.roi_center(event)
        dx, dy = offset
        x0 += dx
        y0 += dy
        half_side = 0.5 * canvas_side
        plt.gca().set_xlim(x0 - half_side, x0 + half_side)
        plt.gca().set_ylim(y0 - half_side, y0 + half_side)
        return collection

    @staticmethod
    def show_display(file_path=None, dpi=100, batch=False):
        """Convenience function to setup the matplotlib canvas for an event display.

        Arguments
        ---------
        file_path : str
            Optional file path to save the image immediately before the plt.show() call.
        """
        plt.gca().set_aspect('equal')
        plt.axis('off')
        if file_path is not None:
            logger.info('Saving event display to %s...', file_path)
            plt.savefig(file_path, dpi=dpi)
        if not batch:
            logger.info('Showing event display, close the window to move to the next one...')
            plt.show()



class xXpolGrid(xHexagonalGrid):

    """XPOL grid.
    """

    def __init__(self, **kwargs):
        """Constructor.
        """
        super().__init__(*XPOL_SIZE, XPOL_PITCH, **kwargs)



class xDisplayArgumentParser(xArgumentParser):

    """Specialized argument parser for the event display and related facilities.

    This is placed here because if needs to be used by both the single-event
    display and the observation carousel.
    """

    def __init__(self, description):
        """Constructor.
        """
        xArgumentParser.__init__(self, description)
        self.add_file()
        self.add_argument('--evtlist', type=str,
            help='path to the auxiliary (Level-2 file) event list')
        self.add_argument('--targetname', type=str, default='N/A',
            help='name of the celestial target')
        self.add_ebounds()
        self.add_seed(default=1)
        self.add_boolean('--clustering', True,
            help='run the DBscan clustering on the events')
        self.add_argument('--numclusters', type=int, default=2,
            help='the number of clusters to be displayed for each event')
        self.add_argument('--clumindensity', type=int, default=5,
            help='the minimum density point for the DBscan clustering')
        self.add_argument('--cluminsize', type=int, default=6,
            help='the minimum cluster size for the DBscan clustering')
        self.add_argument('--resample', type=float, default=None,
            help='the power-law index for resampling events in energy')
        self.add_boolean('--absorption', True,
            help='draw the reconstructed absorption_point')
        self.add_boolean('--barycenter', True,
            help='draw the reconstructed barycenter')
        self.add_boolean('--direction', True,
            help='draw the reconstructed track direction')
        self.add_boolean('--pixpha', False,
            help='indicate the pixel PHA values')
        self.add_boolean('--indices', False,
            help='draw the row and column indices of the readout matrix')
        self.add_argument('--cmap', type=str, default='Reds',
            help='the color map for the pixel values')
        self.add_argument('--cmapoffset', type=int, default=10,
            help='the PHA offset for the color map')
        self.add_argument('--minaxside', type=float, default=2.,
            help='the axis side for the event display')
        self.add_argument('--autostop', type=int, default=None,
            help='stop automatically after a given number of events')
        self.add_boolean('--batch', default=False,
            help='run in batch mode')
        self.add_boolean('--autosave', False,
            help='save the event displays automatically')
        self.add_outfolder(default=IXPEOBSSIM_DATA)
        self.add_argument('--imgformat', type=str, default='png',
            help='the image format for the output files when autosave is True')
        self.add_argument('--imgdpi', type=int, default=250,
            help='resolution of the output image in dot per inches')




def load_event_list(file_path, pivot_energy=8., interactive=False, **kwargs):
    """Load the event data from the Level-2 event list for the purpose of the
    event display---these include, in order: mission elapsed time, energy,
    sky position and Stokes parameters.

    This function has a few other functionalities, other than just loadinf the
    relevant colums from the the Level-2 files, and particularly:

    * resample the input events with a given power-law spectral function;
    * trimming the resampled colums to a target number of events preserving the
      time ordering and covering evengly the entire time span.

    Arguments
    ---------
    file_path : str
        The path to the input Level-2 file.

    pivot_energy : float
        The pivot energy for the resampling of the count spectrum.

    interactive : bool
        If True, show some debug plot with the output (resampled) spectrum.

    kwargs : dict
        The keyword arguments from the xDisplayArgumentParser.
    """
    event_file = xEventFile(file_path)
    logger.info('Total good time in the Level-2 file: %.3f', event_file.total_good_time())
    logger.info('Livetime in the Level-2 file: %.3f', event_file.livetime())
    logger.info('Livetime correction: %.3f', event_file.deadtime_correction())
    resample_index = kwargs.get('resample')
    emin = kwargs.get('emin')
    emax = kwargs.get('emax')
    logger.info('Loading event list from %s...', file_path)
    energy = event_file.energy_data()
    logger.info('Selecting energies in %.2f--%.2f keV...', emin, emax)
    mask = numpy.logical_and(energy >= emin, energy < emax)
    logger.info('Done, %d event(s) out of %d remaining.', mask.sum(), len(mask))
    energy = energy[mask]
    met = event_file.time_data()[mask]
    ra, dec = event_file.sky_position_data()
    ra = ra[mask]
    dec = dec[mask]
    q, u = event_file.stokes_data()
    q = q[mask]
    u = u[mask]
    if resample_index is not None:
        logger.info('Resampling input level-2 data with index %.3f', resample_index)
        mask = numpy.random.uniform(size=len(energy)) <= (energy /  pivot_energy)**resample_index
        logger.info('Done, %d event(s) out of %s remaining.', mask.sum(), len(mask))
        met, energy, ra, dec, q, u = [item[mask] for item in (met, energy, ra, dec, q, u)]
    if interactive:
        # Debug plot for the input energy spectrum.
        plt.figure('Input energy spectrum')
        h = xHistogram1d(numpy.linspace(2., 8., 20)).fill(energy)
        h.plot()
    autostop = kwargs.get('autostop')
    if autostop is not None and autostop < len(met):
        logger.info('Trimming down the L2 columns to the target autostop...')
        # We achieve this by trying and select events uniformly within the range.
        mask = numpy.zeros(len(met), dtype=bool)
        num_events = len(met)
        start_event = num_events // autostop
        idx = numpy.linspace(start_event, num_events - 1, autostop, dtype=int).astype(int)
        mask[idx] = True
        met, energy, ra, dec, q, u = [item[mask] for item in (met, energy, ra, dec, q, u)]
        logger.info('Done, %d event(s) left.', len(met))
    return met, energy, ra, dec, q, u



class xDisplayCard(xTextCard):

    """Specialize text card to display event information.

    The basic idea, here, is that one initializes the card with the EVENTS
    header of a Level-2 file, and then updates the information on an event-by-event
    basis using the set_event_data() hook.
    """

    def __init__(self, target_name, header):
        """Constructor.
        """
        xTextCard.__init__(self)
        self.set_line('Target Name', '%s (obs. %s)' % (target_name, header['OBS_ID']))
        self.set_line('Observation Start', header['DATE-OBS'])
        self.set_line('Observation End', header['DATE-END'])
        self.set_line('Detector Unit', '%s (%s)' % (header['DETNAM'], header['DET_ID']))
        for i in range(5):
            self.set_line('Spacer%d' % i, None)

    def update_cumulative_statistics(self, num_events, emin, emax):
        """Set the card line with the basic cumulative statistics info.
        """
        key = 'Accumulated statistics in %.1f-%.1f keV' % (emin, emax)
        text = '%d events' % num_events
        self.set_line(key, text)

    def set_event_data(self, met, energy, ra, dec, q, u, compact=True):
        """Set the event data.
        """
        self.set_line('Mission elapsed time', met, '%.6f', 's')
        self.set_line('Energy', energy, '%.2f', 'keV')
        if compact:
            self.set_line('Sky position (R. A., Dec.)', '(%.3f, %.3f) decimal degrees' % (ra, dec))
        else:
            self.set_line('Right ascention', ra, '%.3f', 'decimal degrees')
            self.set_line('Declination', dec, '%.3f', 'decimal degrees')
        self.set_line('Stokes parameters (q, u)', '(%.4f, %.4f)' % (q, u))



def display_event(event, grid, threshold, dbscan, file_name=None, padding=False, **kwargs):
    """Single-stop event display.
    """
    draw_kwargs = dict(values=kwargs.get('pixpha'), indices=kwargs.get('indices'),
        min_canvas_side=kwargs.get('minaxside'), zero_sup_threshold=threshold,
        padding=padding, num_clusters=kwargs.get('numclusters'))
    logger.info('Drawing event @ MET %.6f', event.timestamp)
    if kwargs.get('clustering'):
        event.run_clustering(dbscan)
    # Draw the bare event...
    grid.draw_event(event, **draw_kwargs)
    # ... then the reconstruction elements.
    if kwargs.get('absorption'):
        event.recon.draw_absorption_point()
    if kwargs.get('barycenter'):
        event.recon.draw_barycenter()
    if kwargs.get('direction'):
        event.recon.draw_track_direction()
    # Draw the small ixpeobssim signature :-)
    plt.text(0., 0.02, 'Powered by https://github.com/lucabaldini/ixpeobssim',
        transform = plt.gca().transAxes, size='xx-small')
    if kwargs.get('autosave'):
        file_path = os.path.join(kwargs.get('outfolder'), file_name)
    else:
        file_path = None
    grid.show_display(file_path, kwargs.get('imgdpi'), kwargs.get('batch'))

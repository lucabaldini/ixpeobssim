#!/usr/bin/env python
#
# Copyright (C) 2015--2019, the ixpeobssim team.
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

import numbers
import os

import numpy
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.math_ import format_value_error
from ixpeobssim.utils.time_ import met_to_num, met_to_mjd, seconds_to_days


# Hack to support different matplotlib versions.
# MouseButton was introduced in 3.1.1
try:
  from matplotlib.backend_bases import MouseButton
except ImportError:
  class MouseButton:
    LEFT = 1
    MIDDLE = 2
    RIGHT = 3

# A safe place for storing object we do not want to be garbage-collected
STORE_OBJECT_POOL = []

DEFAULT_FIG_WIDTH = 8.

DEFAULT_FIG_HEIGHT = 6.

DEFAULT_FIG_SIZE = (DEFAULT_FIG_WIDTH, DEFAULT_FIG_HEIGHT)

DEFAULT_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
    '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
]

SERIF_FONTS = [
    'DejaVu Serif', 'Bitstream Vera Serif', 'New Century Schoolbook',
    'Century Schoolbook L', 'Utopia', 'ITC Bookman', 'Bookman',
    'Nimbus Roman No9 L', 'Times New Roman', 'Times', 'Palatino', 'Charter',
    'serif'
]

SANS_SERIF_FONTS = [
    'DejaVu Sans', 'Bitstream Vera Sans', 'Lucida Grande', 'Verdana', 'Geneva',
    'Lucid', 'Arial', 'Helvetica', 'Avant Garde', 'sans-serif'
]

CURSIVE_FONTS = [
    'Apple Chancery', 'Textile', 'Zapf Chancery', 'Sand', 'Script MT',
    'Felipa', 'cursive'
]

FANTASY_FONTS = [
    'Comic Sans MS', 'Chicago', 'Charcoal', 'Impact', 'Western', 'Humor Sans',
    'xkcd', 'fantasy'
]

MONOSPACE_FONTS = [
    'DejaVu Sans Mono', 'Bitstream Vera Sans Mono', 'Andale Mono',
    'Nimbus Mono L', 'Courier New', 'Courier', 'Fixed', 'Terminal', 'monospace'
]


# This is a horrible hack to try and support different versions of
# matplotlib---and different versions of cycler.
#
# First try and import cycler.
try:
    from cycler import cycler
    try:
        # Try and use the "new" syntax for cycler 0.10.0 and above.
        COLOR_CYCLER = cycler(color=DEFAULT_COLORS)
    except TypeError:
        # This will give a TypeError for older cycler versions, for which
        # we resort to the old syntax.
        COLOR_CYCLER = cycler('color', DEFAULT_COLORS)
except ImportError:
    logger.warning('Cannot import cycler, falling back to default color wheel...')
    COLOR_CYCLER = None


_COLOR_WHEEL_AR = ['fuchsia', 'green', 'orange', 'dodgerblue', 'magenta']
_COLOR_WHEEL_AR_LEN = len(_COLOR_WHEEL_AR)
def color_wheel_ar(i):
    """Support for a custom color wheel.

    This was originally defined in evt.colorselection as setEnergyColor() and
    was moved here for consistency, see issue #251.
    """
    return _COLOR_WHEEL_AR[i % _COLOR_WHEEL_AR_LEN]


_COLOR_WHEEL_MPR = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
_COLOR_WHEEL_MPR_LEN = len(_COLOR_WHEEL_MPR)
def color_wheel_mpr(i):
    """Support for a custom color wheel.

    This was originally defined in __init__ as xpColor() and
    was moved here for consistency, see issue #252.
    """
    return _COLOR_WHEEL_MPR[i % _COLOR_WHEEL_MPR_LEN]



def du_color(du_id):
    """Return the default color for a give DU.
    """
    return DEFAULT_COLORS[du_id - 1]


def save_gcf(output_folder=IXPEOBSSIM_DATA, file_name=None, file_extensions=('pdf', 'png')):
    """Save the current matplotlib figure.

    If the current figure has a sensible name, it will be used to construct
    the path to the output file---we just make everything lower case and
    replace spaces with `_`.

    Examples
    --------
    >>> from ixpeobssim.utils.matplotlib_ import *
    >>> plt.figure('Test figure')
    >>> # ... do something.
    >>> # This will create `output_folder/test_figure.png`.
    >>> save_gcf('output_folder')

    Arguments
    ---------
    output_folder : string
        The path to the output folder (default to pwd).

    file_name : string
        The figure name (the name of the output files will be name.extension).

    file_extensions : list of strings
        A list of extensions the figure needs to be saved into.

    Returns
    -------
    list
        The list of paths to the file(s) being created by the function.
    """
    file_list = []
    if file_name is None:
        file_name = plt.gcf().get_label().lower().replace(' ', '_')
    assert file_name != ''
    for ext in file_extensions:
        file_path = os.path.join(output_folder, '%s.%s' % (file_name, ext))
        logger.info('Saving current figure to %s...', file_path)
        try:
            plt.savefig(file_path, transparent=False)
            file_list.append(file_path)
        except IOError as e:
            logger.error(e)
    return file_list


def save_all_figures(output_folder, file_extensions=('png', 'pdf')):
    """Save all the figures in memory to a given output folder.
    """
    logger.info('Looping over matplotlib figures and saving all of them...')
    for figid in plt.get_fignums():
        plt.figure(figid)
        save_gcf(output_folder, file_extensions=file_extensions)


def residual_plot(figure_name=None, separation=0.3, padding=0.01,
                  ylabel_offset=-0.085, figsize=DEFAULT_FIG_SIZE):
    """Create a new figure with two axes objects for residual plots.

    Arguments
    ---------
    figure_name : str
        The name of the figure.

    separation : float
        The vertical separation point between the two axes.
    """
    left = rc_param('figure.subplot.left')
    right = rc_param('figure.subplot.right')
    bot = rc_param('figure.subplot.bottom')
    top = rc_param('figure.subplot.top')
    fig = plt.figure(figure_name, figsize=figsize)
    ax1 = fig.add_axes((left, separation + padding, right - left, top - separation))
    ax1.set_xticklabels([])
    ax1.get_yaxis().set_label_coords(ylabel_offset, 0.5)
    ax2 = fig.add_axes((left, bot, right - left, separation - bot - padding))
    plt.sca(ax1)
    ax2.get_yaxis().set_label_coords(ylabel_offset, 0.5)
    return ax1, ax2


def last_line_color(default='black'):
    """Return the color used to draw the last line
    """
    try:
        return plt.gca().get_lines()[-1].get_color()
    except:
        return default


def marker(x, y, **kwargs):
    """Draw a marker at a specified point.
    """
    plt.plot([x], [y], 'o', **kwargs)


def labeled_marker(x, y, label, dx=0, dy=0, **kwargs):
    """Draw a marker and a label at a specified point.
    """
    if kwargs.get('color') is None:
        kwargs['color'] = last_line_color()
    marker(x, y, color=kwargs.get('color'))
    plt.text(x + dx, y + dy, label, **kwargs)


def plot_arrows(grid, field, threshold=0., **kwargs):
    """Plot a 2-dimensional field of arrows.

    This is a lightweight wrapper upon the plt.quiver() method, that is
    taylored to the overlay of arrow fields upon polarizaion degree maps.

    Arguments
    ---------
    grid : array_like
        The underlying grid for the display

    field : array_like or callable
        Anything that can be called on the grid and returns the components of
        the arrow field in the grid points, or an array matching the grid
        that can be used directly

    threshold : float, optional
        Optional threshold on the magnitude of the field at any given point
        (points below the threshold are set to zero)

    kwargs : dict
        All the keyword arguments passed to plt.quiver()
    """
    kwargs.setdefault('color', 'white')
    kwargs.setdefault('alpha', 0.8)
    kwargs.setdefault('transform', plt.gca().get_transform('world'))
    kwargs.setdefault('width', 0.003)
    kwargs.setdefault('headlength', 0.)
    kwargs.setdefault('headwidth', 1.)
    kwargs.setdefault('pivot', 'middle')
    kwargs.setdefault('scale', None)
    xmin, xmax, ymin, ymax = plt.axis()
    x, y = grid
    if callable(field):
        dx, dy = field(x, y)
    else:
        dx, dy = field
    mask = dx**2. + dy**2. > threshold
    if len(mask) > 0:
        plt.gca().quiver(x[mask], y[mask], dx[mask], dy[mask], **kwargs)
    plt.axis([xmin, xmax, ymin, ymax])


def nlog_errorbars(x, y, dy, **kwargs):
    """Plot an errorbar by taking the absolute value of y and flagging the
    negative part with hollow markers.

    This is useful to plot Stokes parameters (that can go negative) in log scale.
    """
    # Create the positive and negative masks for the errorbars.
    _pos = y > 0.
    _neg = numpy.logical_not(_pos)
    y = abs(y)
    # Plot the positive part passing along all the keyword-arguments.
    plt.errorbar(x[_pos], y[_pos], dy[_pos], **kwargs)
    # For the negative part we want to tweak the keyword arguments.
    kwargs.update(dict(color=last_line_color(), markerfacecolor='none'))
    # Note we don't want the label twice.
    kwargs.pop('label', None)
    plt.errorbar(x[_neg], y[_neg], dy[_neg], **kwargs)


def setup_gca(xlabel=None, ylabel=None, xmin=None, xmax=None, ymin=None,
              ymax=None, logx=False, logy=False, grids=False, xticks=None,
              yticks=None, legend=False):
    """Setup the axes for the current plot.
    """
    if logx is True:
        plt.xscale('log')
    if logy is True:
        plt.yscale('log')
    if xticks is not None:
        plt.gca().set_xticks(xticks)
    if yticks is not None:
        plt.gca().set_yticks(yticks)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if xmin is None and xmax is None and ymin is None and ymax is None:
        pass
    else:
        plt.axis([xmin, xmax, ymin, ymax])
    if grids:
        plt.grid(True, which='both')
    if legend:
        plt.legend()


def plot_circle(center, radius, **kwargs):
    """Plot a circle.
    """
    kwargs.setdefault('fill', None)
    circle = plt.Circle(center, radius, **kwargs)
    plt.gca().add_patch(circle)
    return circle


def plot_ellipse(center, width, height, **kwargs):
    """Plot an ellipse.
    """
    kwargs.setdefault('fill', None)
    ellipse = Ellipse(center, width, height, **kwargs)
    plt.gca().add_patch(ellipse)
    return ellipse


def setup_gca_stokes(side=1., pd_grid=numpy.arange(0.1, 1., 0.1),
    pd_grid_label_angle=45., pa_grid_step=30., line_width=0.75, **kwargs):
    """Setup the current axis object in a way that is appropriate for plotting
    Stokes parameters, i.e., for the (Q/I, U/I) phase space.

    Arguments
    ---------
    side : float
        The absolute value of the minimum and maximum Q/I and U/I values for the
        axes ranges.

    pd_grid : array like
        The polarization degrees values for which we plot the reference circles.

    pd_grid_label_angle : float, optional
        The angle (in decimal degrees) at which the pd level labels are rendered.

    pa_grid_step : float
        The azimuthal angle step (in degrees) for the diagonal lines.

    kwargs
    """
    # Set proper defaults for the keyowrd arguments to be passed to the plain setup_gca.
    kwargs.setdefault('xmin', -side)
    kwargs.setdefault('xmax', side)
    kwargs.setdefault('ymin', -side)
    kwargs.setdefault('ymax', side)
    kwargs.setdefault('xlabel', 'Q/I')
    kwargs.setdefault('ylabel', 'U/I')
    setup_gca(**kwargs)
    # Freeze the aspect ratio of the plot.
    plt.gca().set_aspect('equal')
    text_kwargs = dict(ha='center', va='center', backgroundcolor='white', size='small')
    # Draw diagonal lines at the proper angles.
    r = pd_grid[-1]
    for phi in numpy.arange(0., 360., pa_grid_step):
        x = r * numpy.cos(numpy.radians(phi))
        y = r * numpy.sin(numpy.radians(phi))
        plt.plot([0, 2. * x], [0, 2. * y], color='gray', ls='dashed', lw=line_width)
        pa = numpy.degrees(0.5 * numpy.arctan2(y, x))
        # Small hack to move a little bit on the inner side the PA labels
        # at 0, +/- 45 and 90.
        for _pa in (0., 45., -45., 90.):
            if abs(pa - _pa) < 5.:
                x *= 0.9
                y *= 0.9
                break
        plt.text(x, y, '%d$^\\circ$' % numpy.round(pa), **text_kwargs)
    # Draw the contours in polarization degree, i.e., Q^2 + U^2.
    for pd in pd_grid:
        plot_circle((0, 0), pd, color='gray', ls='dashed', lw=line_width)
        if pd_grid_label_angle is not None:
            x = pd * numpy.sin(numpy.radians(pd_grid_label_angle))
            y = pd * numpy.cos(numpy.radians(pd_grid_label_angle))
            plt.text(x, y, '%.2f' % pd, **text_kwargs)


def metplot(met, values, display='dt', *args, **kwargs):
    """
    """
    assert display in ['dt', 'mjd']
    if display == 'dt':
        t = met_to_num(met)
        ylabel = ''
    plt.plot(t, values, *args, **kwargs)
    span = seconds_to_days(met.max() - met.min())
    if span <= 1:
        locator = matplotlib.dates.MinuteLocator()
    elif span <= 31:
        locator = matplotlib.dates.DayLocator()
    else:
        locator = matplotlib.dates.MonthLocator()
    formatter = matplotlib.dates.AutoDateFormatter(locator)
    plt.gca().xaxis.set_major_formatter(formatter)
    setup_gca(xmin=t.min(), xmax=t.max(), ylabel=ylabel)


class xStatBox:

    """Base class describing a text box, to be used for the histogram and fit
    stat boxes.

    In the initial implementation this was wrapped into a simple function,
    but we immediately run into problems in cases where we needed to
    setup a stat box in differen steps before drawing it on the canvas
    (e.g., when a FitModel subclass needs to customize the stat box of the
    base class). The class approach is more flexible, although one needs
    a few more lines of code to add entries and plot the thing.

    Parameters
    ----------
    position : str of tuple
        It can either be a two-element tuple (in which case the argument is
        interpreted as a position in absolute coordinates, with the reference
        corner determined by the alignment flags), or a string in the
        set ['upper left', 'upper right', 'lower left', 'lower rigth'].
        If position is a string, the alignment flags are ignored.

    halign : str
        The horizontal alignment ('left' | 'center' | 'right')

    valign : str
        The vertical alignment ('top' | 'center' | 'bottom')
    """

    HORIZONTAL_PADDING = 0.025
    VERTICAL_PADDING = 0.035
    _left, _right = HORIZONTAL_PADDING, 1 - HORIZONTAL_PADDING
    _bottom, _top = VERTICAL_PADDING, 1 - VERTICAL_PADDING
    POSITION_DICT = {
        'upper left': (_left, _top, 'left', 'top'),
        'upper right': (_right, _top, 'right', 'top'),
        'lower left': (_left, _bottom, 'left', 'bottom'),
        'lower right': (_right, _bottom, 'right', 'bottom')
    }
    DEFAULT_BBOX = dict(boxstyle='round', facecolor='white', alpha=0.75)

    def __init__(self, position='upper left', halign='left', valign='top'):
        """Constructor.
        """
        self.set_position(position, halign, valign)
        self.text = ''

    def set_position(self, position, halign='left', valign='top'):
        """Set the position of the bounding box.
        """
        if isinstance(position, str):
            self.x0, self.y0, self.halign,\
                self.valign = self.POSITION_DICT[position]
        else:
            self.x0, self.y0 = position
            self.halign, self.valign = halign, valign

    def add_entry(self, label, value=None, error=None):
        """Add an entry to the stat box.
        """
        if value is None and error is None:
            self.text += '%s\n' % label
        elif value is not None and error is None:
            try:
                self.text += '%s: %g\n' % (label, value)
            except TypeError:
                self.text += '%s: %s\n' % (label, value)
        elif value is not None and error is not None:
            if error > 0:
                self.text += '%s: %s\n' %\
                    (label, format_value_error(value, error, pm='$\\pm$'))
            else:
                self.text += '%s: %g (frozen)\n' % (label, value)

    def plot(self, **kwargs):
        """Plot the stat box.

        Parameters
        ----------
        **kwargs : dict
            The options to be passed to `plt.text()`
        """
        def set_kwargs_default(key, value):
            """
            """
            if key not in kwargs:
                kwargs[key] = value

        set_kwargs_default('horizontalalignment', self.halign)
        set_kwargs_default('verticalalignment', self.valign)
        set_kwargs_default('bbox', self.DEFAULT_BBOX)
        set_kwargs_default('transform', plt.gca().transAxes)
        plt.text(self.x0, self.y0, self.text.strip('\n'), **kwargs)



class xTextCard(dict):

    """Small class reperesenting a text card.

    This is essentially a dictionary that is capable of plotting itself on
    a matplotlib figure in the form of a multi-line graphic card.
    """

    def __init__(self):
        """Constructor.
        """
        dict.__init__(self)
        self.key_kwargs = dict(color='gray', size='small', ha='left', va='top')
        self.value_kwargs = dict(color='black', ha='left', va='top')

    def set_line(self, key, value, fmt='%s', units=None):
        """Set the value for a given key.

        Arguments
        ---------
        key : str
            The key, i.e., the explanatory text for a given value.

        value : number or str
            The actual value (if None, a blank line will be added).

        fmt : str
            The string format to be used to render the value.

        units : str
            The measurement units for the value.
        """
        self[key] = (value, fmt, units)

    def draw(self, x0=0.1, y0=0.9, line_spacing=0.1, spacing_ratio=0.9):
        """Draw the card.

        Arguments
        ---------
        x0, y0 : float
            The absolute coordinates of the top-left corner of the card.

        line_spacing : float
            The line spacing in units of the total height of the current axes.

        spacing_ratio : float
            The fractional line spacing assigned to the key label.
        """
        key_norm = spacing_ratio / (1. + spacing_ratio)
        value_norm = 1. - key_norm
        for key, (value, fmt, units) in self.items():
            if value is None:
                y0 -= 0.5 * line_spacing
                continue
            plt.text(x0, y0, key, self.key_kwargs)
            y0 -= key_norm * line_spacing
            value = fmt % value
            if units is not None:
                value = '%s %s' % (value, units)
            plt.text(x0, y0, value, self.value_kwargs)
            y0 -= value_norm * line_spacing



class xDraggableColorbar:
    """Interactive colorbar than can be dragged and zoomed.
    Most of the code is taken from:
    http://www.ster.kuleuven.be/~pieterd/python/html/plotting/interactive_colorbar.html

    Input arguments: a standard colorbar and a mappable object (e.g. an image
    created with pyplot.imshow())
    """

    SCROLL_SPEED = 0.05 # Must be (0, 1]
    DRAG_SPEED = 0.05 # Must be (0, 1.].

    def __init__(self, cbar, mappable):
        self.cbar = cbar
        self.mappable = mappable
        self.press = None
        # Create a list of colormaps
        self.cycle = sorted([i for i in dir(plt.cm) if hasattr(getattr(plt.cm,i),'N')])
        # Variable holding the current color map index in the list
        self.index = self.cmap_index()
        # Store the initial configuration for the 'restore' button
        self.initial_vmin = self.cbar.norm.vmin
        self.initial_vmax = self.cbar.norm.vmax
        self.initial_cmap_index = self.index

    @property
    def ax(self):
        return self.cbar.ax

    def set_label(self, zlabel):
        """Set the colorbar label.
        """
        self.cbar.set_label(zlabel)

    def cmap_name(self):
        """ Return the name of the current color map
        """
        return self.mappable.get_cmap().name

    def cmap_index(self):
        """ Return the index of the current color map
        """
        return self.cycle.index(self.mappable.get_cmap().name)

    def show_cmap_name(self):
        """ Show the color map name on top of the color bar
        """
        self.ax.set_title(self.cmap_name())

    def redraw(self):
        """Redraw the colorbar and the related mappable object
        """
        self.cbar.draw_all()
        self.ax.figure.canvas.draw()

    def update_cmap(self, index):
        """Set the color map to the one at the given index
        """
        if self.index == index:
            return
        if index < 0:
            index = len(self.cycle) - 1
        elif index >= len(self.cycle):
            index = 0
        cmap = self.cycle[index]
        self.mappable.set_cmap(cmap)
        self.index = index

    def update_norm(self):
        """ Update the color normalization of the mappable object to reflect
        that of the colorbar
        """
        self.mappable.set_norm(self.cbar.norm)

    def update_limits(self, vmin, vmax):
        """ Update the extrema (vmin and vmax) of the colorbar and update
        the mappable color normalization to it
        """
        self.cbar.norm.vmin = vmin
        self.cbar.norm.vmax = vmax
        self.update_norm()

    def restore(self):
        """Restore the original color map, vmin and vmap
        """
        self.update_cmap(self.initial_cmap_index)
        self.update_limits(self.initial_vmin, self.initial_vmax)

    def connect(self):
        """Connect all the events we need
        """
        self.cidpress = self.ax.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.ax.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.ax.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)
        self.keypress = self.ax.figure.canvas.mpl_connect(
            'key_press_event', self.key_press)

    def on_press(self, event):
        """On button press we will see if the mouse is over the colorbar and
        store the coordinates of the click
        """
        if event.inaxes != self.ax:
            return
        self.press = event.x, event.y

    def key_press(self, event):
        """We connect three butons: 'up' and 'down' to change the colormap and
        'r' to restore the original view (colormap and limits)
        """
        if event.key == 'down':
            self.update_cmap(self.index - 1)
        elif event.key == 'up':
            self.update_cmap(self.index + 1)
        elif event.key == 'r':
            self.restore()
        self.redraw()
        #self.show_cmap_name()

    def on_motion(self, event):
        """On motion we check if the mouse is over the colorbar.
        If the mouse LEFT button is pressed we shift the extrema of the
        colorbar, if the RIGHT button is pressed we zoom in/out
        """
        if self.press is None:
            return
        if event.inaxes != self.ax:
            return
        xprev, yprev = self.press
        dx = event.x - xprev
        dy = event.y - yprev
        self.press = event.x, event.y
        scale = self.cbar.norm.vmax - self.cbar.norm.vmin
        if event.button == MouseButton.LEFT:
            delta = self.DRAG_SPEED * scale * numpy.sign(dy)
            vmin = self.cbar.norm.vmin - delta
            vmax = self.cbar.norm.vmax - delta
        elif (event.button == MouseButton.RIGHT) or \
             (event.button == MouseButton.MIDDLE):
            delta = self.SCROLL_SPEED * scale * numpy.sign(dy)
            vmin = self.cbar.norm.vmin - delta
            vmax = self.cbar.norm.vmax + delta
        else:
            return
        self.update_limits(vmin, vmax)
        self.redraw()

    def on_release(self, event):
        """On release we reset the press data
        """
        self.press = None
        self.update_norm()
        self.redraw()

    def disconnect(self):
        """Disconnect all the stored connection ids
        """
        self.ax.figure.canvas.mpl_disconnect(self.cidpress)
        self.ax.figure.canvas.mpl_disconnect(self.cidrelease)
        self.ax.figure.canvas.mpl_disconnect(self.cidmotion)
        self.ax.figure.canvas.mpl_disconnect(self.keypress)


def draggable_colorbar(mappable=None, cax=None, ax=None, **kwargs):
    """Create a draggable colorbar
    """
    if cax is None and ax is None:
        ax = plt.gca()
    _cbar = plt.colorbar(mappable=mappable, cax=cax, ax=ax, **kwargs)
    cbar = xDraggableColorbar(_cbar, _cbar.mappable)
    cbar.connect()
    # we need to store a reference to the draggable colorbar,
    # otherwise the connection will be destroyed by the garbage collector
    # as soon as the colorbar variable goes out of scope or the canvas
    # loses its focus.
    # See https://matplotlib.org/3.1.0/users/event_handling.html
    STORE_OBJECT_POOL.append(cbar)
    return cbar


def add_slider(val_func, img_obj=None, ax=None, coords=(0.2, 0.05, 0.6, 0.03),
               label=None, valmin=0., valmax=1., num_step=10,
               color='lightgoldenrodyellow', **kwargs):
    """ Add a slider object to a plot. This can be used e.g. to naviagte a
    3-dimensional histogram by plotting a slice of it at a given bin on its
    last axis, which the user can change by scrolling the slider.
    Requires as input a function that returns the new image data based on the
    current value of the slider and the graphical object which is linked
    to the slider (e.g. the one produced by plot() or imgshow()). Default is
    the current image."""
    if ax is None:
        ax = plt.gca()
    if img_obj is None:
        img_obj = plt.gci()
    ax.margins(x=0)
    plt.subplots_adjust(bottom=0.2)
    step = (valmax -  valmin) / num_step
    slider_ax = plt.axes(coords, facecolor=color)
    slider = Slider(slider_ax, label, valmin, valmax, valstep=step, **kwargs)
    def update(val):
        new_data = val_func(slider.val)
        img_obj.set_array(new_data)
        plt.gcf().canvas.draw_idle()

    slider.on_changed(update)
    STORE_OBJECT_POOL.append(slider)




def rc_param(key):
    """Return a given matplotlib configuration property.
    """
    return matplotlib.rcParams[key]


def context_two_by_two(scale=1.9):
    """Setup the current figure for a 2x2 panel.
    """
    _size = (scale*DEFAULT_FIG_WIDTH, scale*DEFAULT_FIG_HEIGHT)
    _rc = {'figure.figsize': _size}
    return matplotlib.rc_context(rc=_rc)


def context_no_grids():
    """Setup the current figure with no grids.
    """
    _rc = {'axes.grid': False}
    return matplotlib.rc_context(rc=_rc)


def _set_rc_param(key, value):
    """Set the value for a single matplotlib parameter.

    The actual command is encapsulated into a try except block because this
    is intended to work across different matplotlib versions. If a setting
    cannot be applied for whatever reason, this will happily move on.
    """
    #logger.debug('Setting %s to %s' % (key, value))
    try:
        matplotlib.rcParams[key] = value
    except Exception:
        pass


def setup():
    """Basic system-wide setup.

    The vast majority of the settings are taken verbatim from the
    matplotlib 2, commit 5285e76:
    https://github.com/matplotlib/matplotlib/blob/master/matplotlibrc.template

    Note that, although this is designed to provide an experience which is as
    consistent as possible across different matplotlib versions, some of the
    functionalities are not implemented in older versions, which is why we wrap
    each parameter setting into a _set_rc_param() function call.
    """
    # http://matplotlib.org/api/artist_api.html#module-matplotlib.lines
    _set_rc_param('lines.linewidth', 1.5)
    _set_rc_param('lines.linestyle', '-')
    _set_rc_param('lines.color', DEFAULT_COLORS[0])
    _set_rc_param('lines.marker', None)
    _set_rc_param('lines.markeredgewidth', 1.0)
    _set_rc_param('lines.markersize', 6)
    _set_rc_param('lines.dash_joinstyle', 'miter')
    _set_rc_param('lines.dash_capstyle', 'butt')
    _set_rc_param('lines.solid_joinstyle', 'miter')
    _set_rc_param('lines.solid_capstyle', 'projecting')
    _set_rc_param('lines.antialiased', True)
    _set_rc_param('lines.dashed_pattern', (2.8, 1.2))
    _set_rc_param('lines.dashdot_pattern', (4.8, 1.2, 0.8, 1.2))
    _set_rc_param('lines.dotted_pattern', (1.1, 1.1))
    _set_rc_param('lines.scale_dashes', True)

    # Markers.
    _set_rc_param('markers.fillstyle', 'full')

    # http://matplotlib.org/api/artist_api.html#module-matplotlib.patches
    _set_rc_param('patch.linewidth', 1)
    _set_rc_param('patch.facecolor', DEFAULT_COLORS[0])
    _set_rc_param('patch.edgecolor', 'black')
    _set_rc_param('patch.force_edgecolor', True)
    _set_rc_param('patch.antialiased', True)

    # Hatches
    _set_rc_param('hatch.color', 'k')
    _set_rc_param('hatch.linewidth', 1.0)

    # Boxplot
    _set_rc_param('boxplot.notch', False)
    _set_rc_param('boxplot.vertical', True)
    _set_rc_param('boxplot.whiskers', 1.5)
    _set_rc_param('boxplot.bootstrap', None)
    _set_rc_param('boxplot.patchartist', False)
    _set_rc_param('boxplot.showmeans', False)
    _set_rc_param('boxplot.showcaps', True)
    _set_rc_param('boxplot.showbox', True)
    _set_rc_param('boxplot.showfliers', True)
    _set_rc_param('boxplot.meanline', False)
    _set_rc_param('boxplot.flierprops.color', 'k')
    _set_rc_param('boxplot.flierprops.marker', 'o')
    _set_rc_param('boxplot.flierprops.markerfacecolor', 'none')
    _set_rc_param('boxplot.flierprops.markeredgecolor', 'k')
    _set_rc_param('boxplot.flierprops.markersize', 6)
    _set_rc_param('boxplot.flierprops.linestyle', 'none')
    _set_rc_param('boxplot.flierprops.linewidth', 1.0)
    _set_rc_param('boxplot.boxprops.color', 'k')
    _set_rc_param('boxplot.boxprops.linewidth', 1.0)
    _set_rc_param('boxplot.boxprops.linestyle', '-')
    _set_rc_param('boxplot.whiskerprops.color', 'k')
    _set_rc_param('boxplot.whiskerprops.linewidth', 1.0)
    _set_rc_param('boxplot.whiskerprops.linestyle', '-')
    _set_rc_param('boxplot.capprops.color', 'k')
    _set_rc_param('boxplot.capprops.linewidth', 1.0)
    _set_rc_param('boxplot.capprops.linestyle', '-')
    _set_rc_param('boxplot.medianprops.color', DEFAULT_COLORS[1])
    _set_rc_param('boxplot.medianprops.linewidth', 1.0)
    _set_rc_param('boxplot.medianprops.linestyle', '-')
    _set_rc_param('boxplot.meanprops.color', DEFAULT_COLORS[2])
    _set_rc_param('boxplot.meanprops.marker', '^')
    _set_rc_param('boxplot.meanprops.markerfacecolor', DEFAULT_COLORS[2])
    _set_rc_param('boxplot.meanprops.markeredgecolor', DEFAULT_COLORS[2])
    _set_rc_param('boxplot.meanprops.markersize', 6)
    _set_rc_param('boxplot.meanprops.linestyle', 'none')
    _set_rc_param('boxplot.meanprops.linewidth', 1.0)

    # http://matplotlib.org/api/font_manager_api.html for more
    _set_rc_param('font.family', 'sans-serif')
    _set_rc_param('font.style', 'normal')
    _set_rc_param('font.variant', 'normal')
    _set_rc_param('font.weight', 'medium')
    _set_rc_param('font.stretch', 'normal')
    _set_rc_param('font.size', 12.0)
    _set_rc_param('font.serif', SERIF_FONTS)
    _set_rc_param('font.sans-serif', SANS_SERIF_FONTS)
    _set_rc_param('font.cursive', CURSIVE_FONTS)
    _set_rc_param('font.fantasy', FANTASY_FONTS)
    _set_rc_param('font.monospace', MONOSPACE_FONTS)

    # http://matplotlib.org/api/artist_api.html#module-matplotlib.text for more
    _set_rc_param('text.color', 'black')

    #http://wiki.scipy.org/Cookbook/Matplotlib/UsingTex
    _set_rc_param('text.usetex', False)
    _set_rc_param('text.hinting', 'auto')
    _set_rc_param('text.hinting_factor', 8)
    _set_rc_param('text.antialiased', True)
    _set_rc_param('mathtext.cal', 'cursive')
    _set_rc_param('mathtext.rm', 'serif')
    _set_rc_param('mathtext.tt', 'monospace')
    _set_rc_param('mathtext.it', 'serif:italic')
    _set_rc_param('mathtext.bf', 'serif:bold')
    _set_rc_param('mathtext.sf', 'sans')
    _set_rc_param('mathtext.fontset', 'stixsans')
    _set_rc_param('mathtext.fallback', 'cm')
    _set_rc_param('mathtext.default', 'it')

    # http://matplotlib.org/api/axes_api.html#module-matplotlib.axes
    _set_rc_param('axes.facecolor', 'white')
    _set_rc_param('axes.edgecolor', 'black')
    _set_rc_param('axes.linewidth', 1.25)
    _set_rc_param('axes.grid', False)
    _set_rc_param('axes.titlesize', 'large')
    _set_rc_param('axes.titlepad', 6.0)
    _set_rc_param('axes.labelsize', 'medium')
    _set_rc_param('axes.labelpad', 4.0)
    _set_rc_param('axes.labelweight', 'normal')
    _set_rc_param('axes.labelcolor', 'black')
    _set_rc_param('axes.axisbelow', 'line')
    _set_rc_param('axes.formatter.limits', (-7, 7))
    _set_rc_param('axes.formatter.use_locale', False)
    _set_rc_param('axes.formatter.use_mathtext', False)
    _set_rc_param('axes.formatter.min_exponent', 0)
    _set_rc_param('axes.formatter.useoffset', True)
    _set_rc_param('axes.formatter.offset_threshold', 4)
    _set_rc_param('axes.spine.left', True)
    _set_rc_param('axes.spine.bottom', True)
    _set_rc_param('axes.spine.top', True)
    _set_rc_param('axes.spine.right', True)
    _set_rc_param('axes.unicode_minus', True)
    _set_rc_param('axes.prop_cycle', COLOR_CYCLER)
    _set_rc_param('axes.autolimit_mode', 'round_numbers')
    _set_rc_param('axes.xmargin', 0.)
    _set_rc_param('axes.ymargin', 0.)
    _set_rc_param('polaraxes.grid', True)
    _set_rc_param('axes3d.grid', True)

    # Dates
    _set_rc_param('date.autoformatter.year', '%Y')
    _set_rc_param('date.autoformatter.month', '%Y-%m')
    _set_rc_param('date.autoformatter.day', '%Y-%m-%d')
    _set_rc_param('date.autoformatter.hour', ' %m-%d %H')
    _set_rc_param('date.autoformatter.minute', '%d %H:%M')
    _set_rc_param('date.autoformatter.second', '%H:%M:%S')
    _set_rc_param('date.autoformatter.microsecond', '%M:%S.%f')

    # see http://matplotlib.org/api/axis_api.html#matplotlib.axis.Tick
    _set_rc_param('xtick.top', False)
    _set_rc_param('xtick.bottom', True)
    _set_rc_param('xtick.major.size', 3.5)
    _set_rc_param('xtick.minor.size', 2)
    _set_rc_param('xtick.major.width', 1.25)
    _set_rc_param('xtick.minor.width', 1.)
    _set_rc_param('xtick.major.pad', 3.5)
    _set_rc_param('xtick.minor.pad', 3.4)
    _set_rc_param('xtick.color', 'k')
    _set_rc_param('xtick.labelsize', 'medium')
    _set_rc_param('xtick.direction', 'out')
    _set_rc_param('xtick.minor.visible', False)
    _set_rc_param('xtick.major.top', True)
    _set_rc_param('xtick.major.bottom', True)
    _set_rc_param('xtick.minor.top', True)
    _set_rc_param('xtick.minor.bottom', True)
    _set_rc_param('ytick.left', True)
    _set_rc_param('ytick.right', False)
    _set_rc_param('ytick.major.size', 3.5)
    _set_rc_param('ytick.minor.size', 2)
    _set_rc_param('ytick.major.width', 1.25)
    _set_rc_param('ytick.minor.width', 1.)
    _set_rc_param('ytick.major.pad', 3.5)
    _set_rc_param('ytick.minor.pad', 3.4)
    _set_rc_param('ytick.color', 'k')
    _set_rc_param('ytick.labelsize', 'medium')
    _set_rc_param('ytick.direction', 'out')
    _set_rc_param('ytick.minor.visible', False)
    _set_rc_param('ytick.major.left', True)
    _set_rc_param('ytick.major.right', True)
    _set_rc_param('ytick.minor.left', True)
    _set_rc_param('ytick.minor.right', True)

    # Grids
    _set_rc_param('grid.color', '#c0c0c0')
    _set_rc_param('grid.linestyle', '--')
    _set_rc_param('grid.linewidth', 0.8)
    _set_rc_param('grid.alpha', 1.0)

    # Legend
    _set_rc_param('legend.loc', 'best')
    _set_rc_param('legend.frameon', True)
    _set_rc_param('legend.framealpha', 0.8)
    _set_rc_param('legend.facecolor', 'inherit')
    _set_rc_param('legend.edgecolor', 0.8)
    _set_rc_param('legend.fancybox', True)
    _set_rc_param('legend.shadow', False)
    _set_rc_param('legend.numpoints', 1)
    _set_rc_param('legend.scatterpoints', 1)
    _set_rc_param('legend.markerscale', 1.0)
    _set_rc_param('legend.fontsize', 'medium')
    _set_rc_param('legend.borderpad', 0.4)
    _set_rc_param('legend.labelspacing', 0.5)
    _set_rc_param('legend.handlelength', 2.0)
    _set_rc_param('legend.handleheight', 0.7)
    _set_rc_param('legend.handletextpad', 0.8)
    _set_rc_param('legend.borderaxespad', 0.5)
    _set_rc_param('legend.columnspacing', 2.0)

    # See http://matplotlib.org/api/figure_api.html#matplotlib.figure.Figure
    _set_rc_param('figure.titlesize', 'large')
    _set_rc_param('figure.titleweight', 'normal')
    _set_rc_param('figure.figsize', DEFAULT_FIG_SIZE)
    _set_rc_param('figure.dpi', 100)
    _set_rc_param('figure.facecolor', 'white')
    _set_rc_param('figure.edgecolor', 'white')
    _set_rc_param('figure.autolayout', False)
    _set_rc_param('figure.max_open_warning', 20)
    _set_rc_param('figure.subplot.left', 0.125)
    _set_rc_param('figure.subplot.right', 0.95)
    _set_rc_param('figure.subplot.bottom', 0.10)
    _set_rc_param('figure.subplot.top', 0.95)
    _set_rc_param('figure.subplot.wspace', 0.2)
    _set_rc_param('figure.subplot.hspace', 0.2)

    # Images
    _set_rc_param('image.aspect', 'equal')
    _set_rc_param('image.interpolation', 'nearest')
    _set_rc_param('image.cmap', 'hot')
    _set_rc_param('image.lut', 256)
    _set_rc_param('image.origin', 'upper')
    _set_rc_param('image.resample', True)
    _set_rc_param('image.composite_image', True)

    # Contour plots
    _set_rc_param('contour.negative_linestyle', 'dashed')
    _set_rc_param('contour.corner_mask', True)

    # Errorbar plots
    _set_rc_param('errorbar.capsize', 0)

    # Histogram plots
    _set_rc_param('hist.bins', 10)

    # Scatter plots
    _set_rc_param('scatter.marker', 'o')

    # Saving figures
    _set_rc_param('path.simplify', True)
    _set_rc_param('path.simplify_threshold', 0.1)
    _set_rc_param('path.snap', True)
    _set_rc_param('path.sketch', None)
    _set_rc_param('savefig.dpi', 'figure')
    _set_rc_param('savefig.facecolor', 'white')
    _set_rc_param('savefig.edgecolor', 'white')
    _set_rc_param('savefig.format', 'png')
    _set_rc_param('savefig.bbox', 'standard')
    _set_rc_param('savefig.pad_inches', 0.1)
    _set_rc_param('savefig.directory', '~')
    _set_rc_param('savefig.transparent', False)

    # Back-ends
    _set_rc_param('pdf.compression', 6)
    _set_rc_param('pdf.fonttype', 3)
    _set_rc_param('svg.image_inline', True)
    _set_rc_param('svg.fonttype', 'path')
    _set_rc_param('svg.hashsalt', None)

    # Key maps
    _set_rc_param('keymap.fullscreen', ('f', 'ctrl+f'))
    _set_rc_param('keymap.home', ('h', 'r', 'home'))
    _set_rc_param('keymap.back', ('left', 'c', 'backspace'))
    _set_rc_param('keymap.forward', ('right', 'v'))
    _set_rc_param('keymap.pan', 'p')
    _set_rc_param('keymap.zoom', 'o')
    _set_rc_param('keymap.save', 's')
    _set_rc_param('keymap.quit', ('ctrl+w', 'cmd+w'))
    _set_rc_param('keymap.grid', 'g')
    _set_rc_param('keymap.grid_minor', 'G')
    _set_rc_param('keymap.yscale', 'l')
    _set_rc_param('keymap.xscale', ('L', 'k'))

    # Animations settings
    _set_rc_param('animation.html', 'none')
    _set_rc_param('animation.writer', 'ffmpeg')
    _set_rc_param('animation.codec', 'h264')
    _set_rc_param('animation.bitrate', -1)
    _set_rc_param('animation.frame_format', 'png')
    _set_rc_param('animation.ffmpeg_path', 'ffmpeg')
    _set_rc_param('animation.ffmpeg_args', '')
    _set_rc_param('animation.mencoder_path', 'mencoder')
    _set_rc_param('animation.mencoder_args', '')
    _set_rc_param('animation.convert_path', 'convert')

    # Miscellanea
    _set_rc_param('timezone', 'UTC')


setup()

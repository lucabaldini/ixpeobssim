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

"""As-run target list.
"""

import datetime

import astropy.coordinates as coord
import astropy.units as u
from matplotlib.dates import date2num, DateFormatter, DayLocator, WeekdayLocator
import numpy

from ixpeobssim.targets.__artl__ import TWG_LIST, COLOR_DICT, _TARGET_DATA,\
    _Y1_ARTL_DATA, SOURCE_LEGEND_HANDLES
from ixpeobssim.utils.matplotlib_ import plt



class xTarget:

    """Small container class representing a celestial target.

    Arguments
    ---------
    name : str
        The source name

    ra : float
        The right ascention of the source in decimal degrees

    dec : float
        The declination of the source in decimal degrees

    twg : TWG entry
        The reference TWG
    """

    def __init__(self, name, ra, dec, twg):
        """Constructor.
        """
        if twg is not None:
            assert twg in TWG_LIST
        self.name = name
        self.ra = ra
        self.dec = dec
        self.twg = twg

    def galactic_coordinates(self):
        """Return the galactic coordinates of the source.
        """
        return coord.SkyCoord(ra=self.ra * u.degree, dec=self.dec * u.degree).galactic

    def _plot(self, x, y, **kwargs):
        """Basic plotting function.
        """
        label = kwargs.pop('label', True)
        ha = kwargs.get('ha', 'center')
        va = kwargs.get('va', 'bottom')
        voff = kwargs.get('voff', 0.02)
        color = COLOR_DICT.get(self.twg)
        rotation = kwargs.get('rotation', 0.)
        plt.scatter(x, y, color=color)
        if label:
            if va == 'bottom':
                y += voff
            else:
                y -= 2. * voff
            plt.text(x, y, self.name, color=color, ha=ha, va=va, rotation=rotation, size='small')

    def plot_equatorial_pos(self, **kwargs):
        """Plot the celestial position.
        """
        ra = coord.Angle(self.ra * u.degree).wrap_at(180. * u.degree).radian
        dec = coord.Angle(self.dec * u.degree).radian
        self._plot(ra, dec, **kwargs)

    def plot_galactic_pos(self, **kwargs):
        """Plot the galactic position.
        """
        coords = self.galactic_coordinates()
        l = -coords.l.wrap_at(180. * u.degree).radian
        b = coords.b.radian
        self._plot(l, b, **kwargs)

    def __str__(self):
        """String formatting.
        """
        return '%s (%s) at (%.5f., %.5f) deg' % (self.name, self.twg, self.ra, self.dec)



class xObservation:

    """Small container representing an IXPE observation.
    """

    _DATETIME_FMT_FULL = '%Y-%m-%dT%H:%M'
    _DATETIME_FMT_TRIM = '%Y-%m-%dT'

    def __init__(self, target_name, exposure, uid, start, stop, note):
        """
        """
        self.target_name = target_name
        self.exposure = self.days_to_ks(exposure)
        self.uid = uid
        self.start_datetime = self._fmt_datetime(start)
        self.end_datetime = self._fmt_datetime(stop)
        self.note = note

    @staticmethod
    def _fmt_datetime(date_str):
        """Format the datetime string into a datetime object.
        """
        try:
            return datetime.datetime.strptime(date_str, xObservation._DATETIME_FMT_FULL)
        except ValueError:
            return datetime.datetime.strptime(date_str, xObservation._DATETIME_FMT_TRIM)

    @staticmethod
    def days_to_ks(days):
        """Convert days to ks.
        """
        return days * 86.400

    def span(self):
        """Return the observation span segment in ks.

        Note this is generally 1.5--2 times longer than the actual exposure, as
        it includes the perdiods in which the target is occulted and or the
        observatory is in the SAA.
        """
        return 1.e-3 * (self.end_datetime - self.start_datetime).total_seconds()

    def __str__(self):
        """String formatting.
        """
        return '%s [%d] %s--%s (%.1f ks)' % (self.target_name, self.uid,
            self.start_datetime, self.end_datetime, self.exposure)



# Static TARGET_DICT object---this can be imported from outside.
TARGET_DICT = {name: xTarget(name, *args) for name, args in _TARGET_DATA.items()}



class xARTL(dict):

    """Small container class representing the as-run target list.
    """

    MAP_FIGSIZE = (14., 9.)

    def __init__(self, data=_Y1_ARTL_DATA):
        """Constructor.
        """
        dict.__init__(self)
        for i, (target_name, *args) in enumerate(data):
            if target_name not in TARGET_DICT:
                msg = 'Target %s not found in _TARGET_DATA, update __artl__.py'
                raise RuntimeError(msg % target_name)
            obs = xObservation(target_name, *args)
            if target_name in self:
                self[target_name].append(obs)
            else:
                self[target_name] = [obs]

    @staticmethod
    def target(target_name):
        """Return the actual target for a given target name.
        """
        return TARGET_DICT[target_name]

    def total_exposure(self, target_name):
        """Return the total exposure for a given target name.
        """
        return sum(obs.exposure for obs in self.get(target_name))

    def _plot_coordinates(self, figname, func, projection='aitoff', **pos_kwargs):
        """
        """
        fig = plt.figure(figname, figsize=self.MAP_FIGSIZE)
        ax = fig.add_subplot(111, projection=projection)
        for target_name in self:
            target = self.target(target_name)
            func(target, **pos_kwargs.get(target_name, {}))
        ax.grid(True)
        plt.xticks(numpy.linspace(-numpy.pi, numpy.pi, 13), labels=[])
        plt.yticks(numpy.linspace(-0.5 * numpy.pi, 0.5 * numpy.pi, 13), labels=[])
        plt.gca().legend(handles=SOURCE_LEGEND_HANDLES, loc=(-0.15, 0.91))

    def plot_equatorial_coodinates(self, projection='aitoff'):
        """Plot the ARTL in equatorial coordinates.
        """
        pos_kwargs = {
            '4U 0142+61': dict(va='top'),
            '4U 1630-472': dict(va='top'),
            '4U 1820-303': dict(ha='left'),
            'Cen X-3': dict(va='top', ha='right'),
            'Cyg X-1': dict(va='top'),
            'Cyg X-2': dict(va='top'),
            'Circinus galaxy': dict(va='top', ha='left'),
            'GS 1826-238': dict(ha='left'),
            'Her X-1': dict(va='top'),
            'NGC 4151': dict(va='top', ha='left'),
            'Sgr A complex': dict(ha='right'),
            'Vela Pulsar': dict(va='top'),
            'XTE J1701-462': dict(ha='left'),
            'PSR B0540-69' : dict(va='bottom'),
            'LMC X-1' : dict(va='top'),
            }
        self._plot_coordinates('IXPE ARTL equatorial', xTarget.plot_equatorial_pos, **pos_kwargs)

    def plot_galactic_coodinates(self, projection='aitoff'):
        """Plot the ARTL in galactic coordinates.
        """
        pos_kwargs = {
            '1RXS J170849.0': dict(va='top', ha='right', rotation=25),
            '4U 0142+61': dict(va='top', ha='right'),
            '4U 1626-67': dict(va='top'),
            '4U 1630-472': dict(ha='left', rotation=25),
            '4U 1820-303': dict(ha='right'),
            'Cas A': dict(va='top'),
            'Circinus galaxy': dict(va='top'),
            'Cen X-3': dict(ha='left'),
            'Cyg X-2': dict(va='top', ha='left'),
            'Cyg X-3': dict(va='top'),
            'GX 301-2': dict(ha='left', rotation=25),
            'GRO J1008-57': dict(va='top'),
            'GRS 1915+105': dict(va='top'),
            'GS 1826-238': dict(va='top', ha='right'),
            'Her X-1': dict(va='top'),
            'Mrk 421': dict(va='top', ha='left'),
            'MSH 15-52': dict(ha='left', rotation=25),
            'NGC 4151': dict(va='top', ha='left'),
            'Vela Pulsar': dict(va='top', ha='left'),
            'Vela X-1': dict(ha='left'),
            'XTE J1701-462': dict(va='top', ha='right', rotation=25),
            'PSR B0540-69' : dict(va='top'),
            'LMC X-1' : dict(va='bottom'),
        }
        self._plot_coordinates('IXPE ARTL galactic', xTarget.plot_galactic_pos, **pos_kwargs)

    def plot(self, y0=0., pad=0., ax=None, label=False, xmin=-numpy.inf, xmax=numpy.inf):
        """Plot the observation timeline.
        """
        if ax is None:
            ax = plt.gca()
        for target_name, obs_list in self.items():
            target = self.target(target_name)
            color = COLOR_DICT[target.twg]
            for obs in obs_list:
                start = date2num(obs.start_datetime)
                if start < xmin or start > xmax:
                    continue
                end = date2num(obs.end_datetime)
                if end < xmin or end > xmax:
                    continue
                ax.hlines(y0, start + pad, end - pad, color=color, lw=25)
                if label:
                    x0 = 0.5 * (start + end)
                    label = '%s' % target_name
                    #if obs.segment_number > 0:
                    #    label += ' (s%d)' % obs.segment_number
                    ax.text(x0, y0 - 0.0275, label, color='black', ha='left', va='bottom', size='small',
                        rotation=45.)
                elif obs.span() > 1300:
                    x0 = 0.5 * (start + end)
                    ax.text(x0, y0, '%s' % target_name, color='white', ha='center', va='center', size='small')

    def plot_inset(self, start_date, end_date, x, y, dy=0.15, xscale=0.0061):
        """Plot a timeline inset.
        """
        _start = date2num(numpy.datetime64(start_date))
        _end = date2num(numpy.datetime64(end_date))
        dx = (_end - _start) * xscale
        axins = plt.gca().inset_axes([x, y, dx, dy])
        axins.set_xlim(_start, _end)
        axins.set_ylim(0.1, -0.1)
        axins.spines['top'].set_visible(False)
        if y > 0.5:
            axins.spines['right'].set_visible(False)
        if y > 0.75:
            axins.spines['left'].set_visible(False)
        axins.patch.set_facecolor('none')
        axins.xaxis_date()
        axins.yaxis.set_ticks([])
        axins.xaxis.set_minor_locator(DayLocator())
        axins.xaxis.set_major_locator(WeekdayLocator())
        axins.xaxis.set_major_formatter(DateFormatter('%b %d'))
        plt.gca().indicate_inset_zoom(axins, edgecolor='black')
        self.plot(y0=0.06, pad=0.05, ax=axins, label=True, xmin=_start, xmax=_end)

    def plot_timeline(self):
        """
        """
        plt.figure('IXPE as-run target list', figsize=(18., 10.))
        self.plot()
        plt.gca().xaxis_date()
        plt.gca().yaxis.set_ticks([])
        plt.axis([None, None, -1., 1.])
        plt.tight_layout()
        self.plot_inset('2022-01-29 12:35', '2022-03-23 01:43', 0.0075, 0.85)
        self.plot_inset('2022-03-24 01:58', '2022-05-14 12:27', 0.0075, 0.25)
        self.plot_inset('2022-05-14 12:28', '2022-06-21 21:04', 0.31, 0.65)
        self.plot_inset('2022-07-07 00:00', '2022-09-19 05:00', 0.32, 0.05)
        self.plot_inset('2022-09-19 05:00', '2022-11-27T13:10', 0.54, 0.85)
        self.plot_inset('2022-11-27T13:10', '2022-12-29T17:45', 0.78, 0.25)
        plt.gca().legend(handles=SOURCE_LEGEND_HANDLES, loc=(0.05, 0.025))

    def __str__(self):
        """String formatting.
        """
        text = ''
        for target_name, obs_list in self.items():
            text += '* %s: %.1f ks total\n' % (target_name, self.total_exposure(target_name))
            for obs in obs_list:
                text += '  - %s\n' % obs
        return text



if __name__ == '__main__':
    artl = xARTL()
    print(artl)
    artl.plot_equatorial_coodinates()
    artl.plot_galactic_coodinates()
    artl.plot_timeline()
    plt.show()

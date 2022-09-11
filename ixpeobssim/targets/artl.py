#!/urs/bin/env python
#
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
import os

import numpy

from ixpeobssim import IXPEOBSSIM_TARGETS_DATA, IXPEOBSSIM_DATA
from ixpeobssim.targets.ltp import xObservationSegment, TARGET_DICT, COLOR_DICT,\
    SOURCE_LEGEND_HANDLES
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, save_all_figures

from matplotlib.dates import date2num, DateFormatter, DayLocator, WeekdayLocator


#
# Name	Exposure [d]	Start	Stop
#
_ARTL_DATA = (
    ('Cas A', 11.57, '2022-01-11T11:45', '2022-01-29T12:34', ''),
    ('Cen X-3',	1.16, '2022-01-29T13:10', '2022-01-31T07:10', ''),
    ('4U 0142+61',	6.94, '2022-01-31T07:36', '2022-02-14T23:48', '1 of 2'),
    ('Cen A', 1.16, '2022-02-15T00:38', '2022-02-17T13:33', ''),
    ('Her X-1', 2.31, '2022-02-17T14:11', '2022-02-21T15:50', '1 of 2'),
    ('Crab', 0.58, '2022-02-21T16:34', '2022-02-22T18:23', '1 of 2'),
    ('Her X-1', 1.16, '2022-02-22T19:07', '2022-02-24T19:24', '2 of 2'),
    ('4U 0142+61',	1.97, '2022-02-24T19:47', '2022-02-27T18:51', '2 of 2'),
    ('Sgr A complex', 4.05, '2022-02-27T19:37', '2022-03-06T23:47', '1 of 2'),
    ('Crab', 0.58, '2022-03-07T00:38', '2022-03-08T02:17', '2 of 2'),
    ('Mrk 501', 1.16, '2022-03-08T02:59', '2022-03-10T08:30', ''),
    ('Sgr A complex', 7.52, '2022-03-10T08:30', '2022-03-23T01:42', '2 of 2'),
    ('4U 1626-67', 2.31, '2022-03-24T01:59', '2022-03-27T05:19', ''),
    ('Mrk 501', 1.16, '2022-03-27T05:59', '2022-03-29T07:03', ''),
    ('GS 1826-238',	1.16, '2022-03-29T07:25', '2022-03-31T08:57', ''),
    ('S5 0716+714',	4.63, '2022-03-31T09:42', '2022-04-05T19:29', ''),
    ('Vela Pulsar', 5.79, '2022-04-05T20:11', '2022-04-15T18:02', '1 of 2'),
    ('Vela X-1', 3.47, '2022-04-15T18:11', '2022-04-21T12:17', ''),
    ('Vela Pulsar', 5.24, '2022-04-21T12:24', '2022-04-30T10:02', '2 of 2'),
    ('Cyg X-2',	1.16, '2022-04-30T11:04', '2022-05-02T11:07', ''),
    ('Cyg X-2', 0.57, '2022-05-02T11:11', '2022-05-03T11:13', 'off-set'),
    ('1ES 1959+650', 0.58,	'2022-05-03T11:29', '2022-05-04T09:49', ''),
    ('Mrk 421',	1.16, '2022-05-04T10:11', '2022-05-06T10:50', ''),
    ('BL Lac', 4.63, '2022-05-06T11:30', '2022-05-14T12:27', ''),
    ('MCG-5-23-16', 0.58, '2022-05-14T13:17', '2022-05-15T14:55', '1 of 2'),
    ('Cyg X-1',	3.47, '2022-05-15T15:45', '2022-05-21T17:52', ''),
    ('MCG-5-23-16', 5.21, '2022-05-21T18:42', '2022-05-31T03:47', '2 of 2'),
    ('3C 454.3', 1.16, '2022-05-31T04:36', '2022-06-02T08:40', ''),
    ('3C 273', 1.16, '2022-06-02T08:53', '2022-06-04T10:48', ''),
    ('Mrk 421',	1.16, '2022-06-04T11:05', '2022-06-06T10:57', ''),
    ('1ES 1959+650', 0.58, '2022-06-06T11:20', '2022-06-07T08:38', 'adjust'),
    ('Mrk 421',	1.16, '2022-06-07T09:01', '2022-06-09T09:39', ''),
    ('1ES 1959+650', 2.46, '2022-06-09T10:02', '2022-06-12T20:25', ''),
    ('3C 279', 3.24, '2022-06-12T21:05', '2022-06-18T20:21', ''),
    ('Cyg X-1', 1.16, '2022-06-18T21:03', '2022-06-21T21:03', 'ToO'),
    ('Tycho', 11.57, '2022-06-21T21:22', '2022-07-07T00:00', ''),
    ('Cen X-3', 2.31, '2022-07-04T06:22', '2022-07-07T12:44', 'ToO'),
    ('BL Lac', 1.44, '2022-07-07T13:35', '2022-07-09T23:12', ''),
    ('Mrk 501', 1.16, '2022-07-09T23:33', '2022-07-12T00:40', ''),
    ('Circinus galaxy', 9.26, '2022-07-12T01:23', '2022-07-25T02:05', ''),
    ('SN 1006', 2.31, '2022-07-25T02:20', '2022-07-29T02:19', 'Seg 1/2 Wrong coordinates'),
    ('GX 301-2', 3.47, '2022-07-29T02:35', '2022-08-03T06:01', ''),
    ('SN 1006', 1.27, '2022-08-03T06:16', '2022-08-05T14:18', 'Seg 2/2: Obs UID=1001599 Wrong coordinates'),
    ('SN 1006', 7.99, '2022-08-05T14:21', '2022-08-19T12:32', ''),
    ('X Persei', 1.61, '2022-08-19T13:22', '2022-08-22T05:39', 'Seg 1/2 interrupted by ToO'),
    ('4U 1630-472', 6.80, '2022-08-22T06:28', '2022-08-23T22:37', 'ToO')
    )





class xTargetList(dict):

    """Small container class representing the as-run target list as
    distributed through https://ixpe.msfc.nasa.gov/for_scientists/asrun.html
    """

    _DATETIME_FMT = '%Y-%m-%dT%H:%M'

    def __init__(self):
        """
        """
        dict.__init__(self)
        for i, (target_name, _, start_date, end_date, _) in enumerate(_ARTL_DATA):
            start_date = datetime.datetime.strptime(start_date, self._DATETIME_FMT)
            end_date = datetime.datetime.strptime(end_date, self._DATETIME_FMT)
            obs = xObservationSegment(target_name, start_date, end_date, 0)
            if target_name in self:
                self[target_name].append(obs)
            else:
                self[target_name] = [obs]

    def total_observing_time(self, target_name):
        """
        """
        return sum(obs.duration() for obs in self.get(target_name))

    def plot(self, y0=0., pad=0., ax=None, label=False, xmin=-numpy.inf, xmax=numpy.inf):
        """
        """
        if ax is None:
            ax = plt.gca()
        for target_name, obs_list in self.items():
            target = TARGET_DICT[target_name]
            color = COLOR_DICT[target.twg]
            for obs in obs_list:
                start = date2num(obs.start_date)
                if start < xmin or start > xmax:
                    continue
                end = date2num(obs.end_date)
                if end < xmin or end > xmax:
                    continue
                ax.hlines(y0, start + pad, end - pad, color=color, lw=25)
                if label:
                    x0 = 0.5 * (start + end)
                    label = '%s' % target_name
                    if obs.segment_number > 0:
                        label += ' (s%d)' % obs.segment_number
                    ax.text(x0, y0 - 0.0275, label, color='black', ha='left', va='bottom', size='small',
                        rotation=45.)
                elif obs.duration() > 1300:
                    x0 = 0.5 * (start + end)
                    ax.text(x0, y0, '%s' % target_name, color='white', ha='center', va='center', size='small')


    def __str__(self):
        """String formatting.
        """
        text = ''
        for target_name, obs_list in self.items():
            text += '* %s: %.1f ks total\n' % (target_name, self.total_observing_time(target_name))
            for obs in obs_list:
                text += '  - %s\n' % obs
        return text



ARTL = xTargetList()


def _plot_timeline_inset(start_date, end_date, x, y, dy=0.15, target_list=ARTL):
    """
    """
    _start = date2num(numpy.datetime64(start_date))
    _end = date2num(numpy.datetime64(end_date))
    dx = (_end - _start) * 0.010
    axins = plt.gca().inset_axes([x, y, dx, dy])
    axins.set_xlim(_start, _end)
    axins.set_ylim(0.1, -0.1)
    axins.spines['top'].set_visible(False)
    if x > 0.8:
        axins.spines['left'].set_visible(False)
    elif y > 0.5:
        axins.spines['right'].set_visible(False)
        axins.spines['left'].set_visible(False)
    axins.patch.set_facecolor('none')
    axins.xaxis_date()
    axins.yaxis.set_ticks([])
    axins.xaxis.set_minor_locator(DayLocator())
    axins.xaxis.set_major_locator(WeekdayLocator())
    axins.xaxis.set_major_formatter(DateFormatter('%b %d'))
    plt.gca().indicate_inset_zoom(axins, edgecolor='black')
    target_list.plot(y0=0.06, pad=0.05, ax=axins, label=True, xmin=_start, xmax=_end)


def plot_timeline():
    """
    """
    plt.figure('IXPE as-run target list', figsize=(18., 10.))
    ARTL.plot()
    plt.gca().xaxis_date()
    plt.gca().yaxis.set_ticks([])
    plt.axis([None, None, -1., 1.])
    plt.tight_layout()
    _plot_timeline_inset('2022-01-29 12:35', '2022-03-23 01:43', 0.0075, 0.82)
    _plot_timeline_inset('2022-03-24 01:58', '2022-05-14 12:27', 0.0075, 0.25)
    _plot_timeline_inset('2022-05-14 12:28', '2022-06-21 21:04', 0.55, 0.65)
    _plot_timeline_inset('2022-07-07 00:00', '2022-08-24 00:00', 0.50, 0.05)
    plt.gca().legend(handles=SOURCE_LEGEND_HANDLES, loc=(0.05, 0.025))



if __name__ == '__main__':
    print(ARTL)
    plot_timeline()
    save_all_figures(IXPEOBSSIM_DATA)
    plt.show()

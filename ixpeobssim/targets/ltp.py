#!/urs/bin/env python
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

"""Baseline long-term observation plan.
"""

import os
import datetime
from enum import Enum, unique

import astropy.coordinates as coord
import astropy.units as u
import numpy
from matplotlib.lines import Line2D
from matplotlib.dates import date2num, DateFormatter, DayLocator, WeekdayLocator

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim import IXPEOBSSIM_TARGETS_DATA
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, save_all_figures




@unique
class TWG(Enum):

    """Definition of the TWGs.
    """

    PWN = 'PWN and radio pulsars'
    SNR = 'SNR'
    ASMBH = 'Accreting stellar-mass BH'
    AWDNS = 'Accreting WD and NS'
    MAGNETAR = 'Magnetars'
    RQAGN = 'Radio-quiet AGN and Sgr A*'
    BLAZAR = 'Blazars and radiogalaxies'


TWGS = list(TWG)


_COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']
COLOR_DICT = {cls: _COLORS[i] for i, cls in enumerate(TWGS)}

SOURCE_LEGEND_HANDLES = []
for i, cls in enumerate(TWG):
    handle = Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_DICT[cls],
        markersize=10., label=cls.value)
    SOURCE_LEGEND_HANDLES.append(handle)



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

    min_flux : float
        The minimum integral X-ray flux between 2 and 10 keV in cgs

    min_flux : float, optional
        The maximum integral X-ray flux between 2 and 10 keV in cgs
    """

    def __init__(self, name, ra, dec, twg, visibility_frac, min_flux, max_flux=None):
        """Constructor.
        """
        assert twg in TWGS
        self.name = name
        self.ra = ra
        self.dec = dec
        self.twg = twg
        self.visibility_frac = visibility_frac
        self.min_flux = min_flux
        if max_flux is None:
            self.max_flux = self.min_flux
        else:
            self.max_flux = max_flux

    def average_flux(self):
        """Return the geometric mean of the minimum and maximum fluxes.
        """
        return (self.min_flux * self.max_flux)**0.5

    def plot_celestial_pos(self, **kwargs):
        """Plot the celestial position.
        """
        label = kwargs.pop('label', True)
        ha = kwargs.get('ha', 'center')
        va = kwargs.get('va', 'bottom')
        voff = kwargs.get('voff', 0.02)
        color = COLOR_DICT.get(self.twg)
        ra = coord.Angle(self.ra * u.degree).wrap_at(180. * u.degree).radian
        dec = coord.Angle(self.dec * u.degree).radian
        plt.scatter(ra, dec, color=color)
        if label:
            if va == 'bottom':
                dec += voff
            else:
                dec -= 2. * voff
            plt.text(ra, dec, self.name, color=color, ha=ha, va=va, size='small')

    def plot_galactic_pos(self, **kwargs):
        """Plot the galactic position.
        """
        label = kwargs.pop('label', True)
        ha = kwargs.get('ha', 'center')
        va = kwargs.get('va', 'bottom')
        voff = kwargs.get('voff', 0.02)
        color = COLOR_DICT.get(self.twg)
        gc = coord.SkyCoord(ra=self.ra * u.degree, dec=self.dec * u.degree).galactic
        l = -gc.l.wrap_at(180. * u.degree).radian
        b = gc.b.radian
        plt.scatter(l, b, color=color)
        if label:
            if va == 'bottom':
                b += voff
            else:
                b -= 2. * voff
            plt.text(l, b, self.name, color=color, ha=ha, va=va, size='small')


    def plot(self, livetime, flux_range=True, **kwargs):
        """Plot the target.
        """
        label = kwargs.pop('label', True)
        ha = kwargs.get('ha', 'center')
        va = kwargs.get('va', 'bottom')
        voff = kwargs.get('voff', 1.03)
        color = COLOR_DICT.get(self.twg)
        x0, y0 = self.average_flux(), livetime / 86.4
        plt.plot(x0, y0, 'o', color=color, markersize=7.5)
        if flux_range and self.min_flux != self.max_flux:
            plt.hlines(y0, self.min_flux, self.max_flux, color=color, lw=1., ls='solid')
        if label:
            if va == 'bottom':
                y0 *= voff
            else:
                y0 /= voff**2.
            plt.text(x0, y0, self.name, color=color, ha=ha, va=va, size='small')



class xTargetDict(dict):

    """Small container class representing the targets in the long-term plan.
    """

    def add(self, *args, **kwargs):
        """
        """
        target_name = args[0]
        self[target_name] = xTarget(*args, **kwargs)



#
# Target dictionary.
#
# The target names have been taken from the long-term plan; the source coordinates
# have been retrieved via xpsrcoords.py, and the average annual visibity has been
# calculated with xpvisibility.py.
#
# The source class and the X-ray flux are from the spreadsheet that Roger Romani
# sent me on December 20, 2021. (The comment at the end of the line is the
# net observation time indicated in the spreadsheet.)
#
TARGET_DICT = xTargetDict()
TARGET_DICT.add('Cas A', 350.85, 58.815, TWG.SNR, .629, 1.5e-9) #1000 ks
TARGET_DICT.add('Cen A', 201.36506288, -43.01911267, TWG.BLAZAR, .566, 1.e-9) #200 ks
TARGET_DICT.add('Her X-1', 254.45754617, 35.34235762, TWG.AWDNS, .551, 0.1e-10, 50.e-10) #400 ks
TARGET_DICT.add('4U 0142+61', 26.59253, 61.75106, TWG.MAGNETAR, .652, 67.9e-12) #1000 ks
TARGET_DICT.add('Crab', 83.63308333, 22.0145, TWG.PWN, .536, 2.5e-8) #100 ks
TARGET_DICT.add('Sgr A complex', 266.41681662, -29.00782497, TWG.RQAGN, .543, 2.e-11) #1000 ks
TARGET_DICT.add('Mrk 501', 253.46756952, 39.76016915, TWG.BLAZAR, .559, 0.5e-11, 8.e-11) #300 ks
TARGET_DICT.add('X Persei', 58.84615783, 31.04584604, TWG.AWDNS, .545, 1.e-10, 5.e-10) #250 ks
TARGET_DICT.add('GS 1826-238', 277.3675, -23.79694444, TWG.AWDNS, .538, 1.e-9, 6.e-9) #100 ks
TARGET_DICT.add('S5 0716+714', 110.47270204, 71.34343391, TWG.BLAZAR, .869, 1.e-11) #400 ks
TARGET_DICT.add('GRS 1915+105', 288.798149, 10.945807, TWG.ASMBH, .530, 1.8e-9, 4.8e-8) #250 ks
TARGET_DICT.add('Vela Pulsar', 128.5, -45.83333333, TWG.PWN, .573, 7.e-12) #1000 ks
TARGET_DICT.add('Mrk 421', 166.113808, 38.20883287, TWG.BLAZAR, .556, 2.e-11, 500.e-11) #300 ks
TARGET_DICT.add('Cyg X-1', 299.59031591, 35.20160625, TWG.ASMBH, .551, 2.4e-11, 2.9e-8) #300 ks
TARGET_DICT.add('MCG-5-23-16', 146.917319, -30.948734, TWG.RQAGN, .545, 0.8e-10, 1.e-10) #500 ks
TARGET_DICT.add('Vela X-1', 135.52858781, -40.55469345, TWG.AWDNS, .560, 0.1e-10, 200.e-10) #300 ks
TARGET_DICT.add('BL Lac', 330.68038074, 42.27777178, TWG.BLAZAR, .564, 1.e-11, 8.e-11) #400 ks
TARGET_DICT.add('3C 454.3', 343.49061658, 16.14821142, TWG.BLAZAR, .533, 0.5e-11, 10.e-11) #200 ks
TARGET_DICT.add('3C 273', 187.27791535, 2.05238857, TWG.BLAZAR, .529, 5.e-11, 40.e-11) #200 ks
TARGET_DICT.add('Cyg X-2', 326.17147688, 38.32140718, TWG.AWDNS, .556, 0.5e-8, 2.e-8) #100 ks
TARGET_DICT.add('1ES 1959+650', 299.99938562, 65.14851421, TWG.BLAZAR, .689, 1.e-11, 150.e-11) #150 ks
TARGET_DICT.add('3C 279', 194.04652688, -5.7893144, TWG.BLAZAR, .529, 1.e-11, 10.e-11) #200 ks
TARGET_DICT.add('Tycho', 6.308299, 64.144305, TWG.SNR, .677, 2.1e-10) #1000 ks
TARGET_DICT.add('Mrk 501', 253.46756952, 39.76016915, TWG.BLAZAR, .559, 0.5e-11, 8.e-11) #300 ks
TARGET_DICT.add('IC 4329A', 207.33025303, -30.3095067, TWG.RQAGN, .544, 8.5842e-11) #500 ks
TARGET_DICT.add('1ES 0229+200', 38.20256456, 20.2881926, TWG.BLAZAR, .535, 1.e-11, 5.e-11) #400 ks
TARGET_DICT.add('J0211+1051', 32.804905, 10.859666, TWG.BLAZAR, .530, 1.e-11) #400 ks
TARGET_DICT.add('GX 301-2', 186.65650295, -62.77034989, TWG.AWDNS, .662, 0.1e-10, 100.e-10) #300 ks
TARGET_DICT.add('Circinus galaxy', 213.291458, -65.339222, TWG.RQAGN, .692, 1.4e-11) #800 ks
TARGET_DICT.add('MSH15-52', 228.32083333, -59.08166667, TWG.PWN, .630, 8.e-12) #1500 ks
TARGET_DICT.add('1RXS J170849.0', 257.20416667, -40.15277778, TWG.MAGNETAR, .560, 24.3e-12) #1000 ks
TARGET_DICT.add('GX9+9', 262.93404167, -16.96136111, TWG.AWDNS, .533, 0.3e-8, 26.5e-8) #100 ks
TARGET_DICT.add('4U 1626-67', 248.06996, -67.46091, TWG.AWDNS, .727, 2.e-10, 5.e-10) #200 ks
TARGET_DICT.add('Cen X-3', 170.31288349, -60.62378542, TWG.AWDNS, .642, 5.e-10)
TARGET_DICT.add('SN 1006', 225.59208, -42.09694, TWG.SNR, .648, None) #1000 ks
TARGET_DICT.add('4U 1630-472', 248.506708, -47.393, TWG.ASMBH, 0.664, None)



class xObservationSegment:

    """Class describing an observation segment.
    """

    def __init__(self, target_name, start_date, end_date, segment_number):
        """Constructor.
        """
        self.target_name = target_name
        self.start_date = start_date
        self.end_date = end_date
        self.segment_number = segment_number

    def duration(self):
        """Return the duration of the observation segment in ks.
        """
        return 1.e-3 * (self.end_date - self.start_date).total_seconds()

    def __str__(self):
        """String formatting.
        """
        return '%s [%d] %s--%s (%.1f ks)' % (self.target_name, self.segment_number,
            self.start_date, self.end_date, self.duration())



class xLongTermPlan(dict):

    """Small container class representing the long-term observing plan as
    distributed through https://ixpe.msfc.nasa.gov/for_scientists/TARGET_DICT.html
    """

    def __init__(self):
        """
        """
        dict.__init__(self)
        ltp_data, self.last_update = self._parse_text_file()
        for i, (target_name, start_date, segment_number) in enumerate(ltp_data):
            try:
                _, end_date, _ = ltp_data[i + 1]
            except IndexError:
                end_date = start_date + datetime.timedelta(days=1)
            obs = xObservationSegment(target_name, start_date, end_date, segment_number)
            if target_name in self:
                self[target_name].append(obs)
            else:
                self[target_name] = [obs]

    def total_observing_time(self, target_name):
        """
        """
        return sum(obs.duration() for obs in self.get(target_name))

    @staticmethod
    def _parse_text_file():
        """Parse a text dump of the long-term plan and extract the information
        in a usable form.
        """
        file_path = os.path.join(IXPEOBSSIM_TARGETS_DATA, 'ixpe_ltp.txt')
        data = []
        last_update = None
        with open(file_path) as input_file:
            accumulating = False
            for line in input_file:
                line = line.strip()
                if line == 'SG	Name	Start	buffer':
                    line = input_file.readline().strip()
                    accumulating = True
                elif line == '':
                    accumulating = False
                if accumulating:
                    segment, *pieces, start_date, _ = line.split()
                    segment = int(segment)
                    start_date = datetime.datetime.strptime(start_date, '%Y-%m-%d')
                    target_name = ' '.join(pieces)
                    data.append((target_name, start_date, segment))
                if 'Author/Curator' in line:
                    last_update = line.split('|')[0]
        if last_update is None:
            logger.error('Could not parse last update.')
        return data, last_update

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



LTP = xLongTermPlan()

def plot_last_update(x=0.99, y=0.05):
    """Plot a last update text label.
    """
    plt.text(x, y, 'LTP updated on %s' % LTP.last_update, ha='right',
        transform=plt.gca().transAxes)


def plot_mdp(flux_range=True):
    """
    """
    plot_kwargs = {
        'Circinus galaxy': dict(va='top'),
        '1RXS J170849.0': dict(va='top'),
        'IC 4329A': dict(va='top'),
        'GRS 1915+105': dict(va='top'),
        '3C 279': dict(va='top'),
        'Cyg X-1': dict(va='top'),
        'GX 301-2': dict(va='top'),
        'Vela X-1': dict(ha='left'),
        '1ES 0229+200': dict(va='top'),
        'MSH15-52': dict(va='top')
    }
    if flux_range:
        plt.figure('IXPE LTP MDP', figsize=(12., 8.))
    else:
        plt.figure('IXPE LTP MDP simple', figsize=(12., 8.))
    for target_name, obs_list in LTP.items():
        target = TARGET_DICT.get(target_name)
        if target is None:
            print('---> Cannot find %s...' % target_name)
            continue
        obs_time = LTP.total_observing_time(target_name)
        obs_livetime = obs_time * target.visibility_frac
        print('%s: %.1f ks wall-time (%.1f ks livetime)' % (target_name, obs_time, obs_livetime))
        target.plot(obs_livetime, label=True, flux_range=flux_range, **plot_kwargs.get(target_name, {}))
    # Our level-1 requirement on the polarization sensitivity reads
    # “MDP (99% confidence) <= 5.5 % for a point source with E**-2 photon spectrum
    # and a 2-8 keV flux of 10--11 ergs/cm2 sec for a 10 day integration.”
    x = numpy.geomspace(1.e-12, 1.e-7, 100)
    for mdp in (30., 10., 3., 1., 0.3):
        y = (5.5 / mdp) * 10. * (x / 1.e-11)**-0.5
        plt.plot(x, y, color='black', ls='dashed', zorder=0)
        if mdp == 30:
            y0 = 3.
        elif mdp == 10:
            y0 = 8.
        else:
            y0 = 15.
        x0 = 1.2 * 1.e-11 * ((5.5 / mdp) * 10. / y0)**2.
        plt.text(x0, y0, 'MDP$_{99}$ = %.1f%%' % mdp, rotation=-45., ha='right')
    setup_gca(xmin=1.e-12, xmax=1.e-7, ymin=0.5, ymax=30., logx=True, logy=True,
        grids=True, ylabel='Effective observation time [days]',
        xlabel='$F_{2\\mathrm{-}10~\\mathrm{keV}}$ [erg cm$^{-1}$ s$^{-1}$]')
    plt.legend(handles=SOURCE_LEGEND_HANDLES)
    plot_last_update()


def plot_celestial(projection='aitoff'):
    """
    """
    plot_kwargs = {
        'Circinus galaxy': dict(va='top'),
        'Vela Pulsar': dict(va='top'),
        'Cyg X-2': dict(va='top'),
        'Her X-1': dict(va='top'),
        '4U 0142+61': dict(va='top')
    }
    fig = plt.figure('IXPE LTP celestial', figsize=(12., 8.))
    ax = fig.add_subplot(111, projection=projection)
    for target in TARGET_DICT.values():
        target.plot_celestial_pos(**plot_kwargs.get(target.name, {}))
    ax.grid(True)
    plt.gca().legend(handles=SOURCE_LEGEND_HANDLES, loc=(-0.15, 0.91))
    plot_last_update(y=-0.05)


def plot_galactic(projection='aitoff'):
    """
    """
    plot_kwargs = {
        'Mrk 421': dict(va='top'),
        'Her X-1': dict(va='top'),
        '4U 0142+61': dict(va='top'),
        'Cas A': dict(va='top'),
        'Cyg X-2': dict(va='top'),
        'GRS 1915+105': dict(va='top'),
        'GS 1826-238': dict(va='top'),
        '4U 1626-67': dict(va='top'),
        'Circinus galaxy': dict(va='top'),
        'Vela Pulsar': dict(va='top'),
        '1RXS J170849.0': dict(va='top'),
    }
    fig = plt.figure('IXPE LTP galactic', figsize=(12., 8.))
    ax = fig.add_subplot(111, projection=projection)
    for target in TARGET_DICT.values():
        target.plot_galactic_pos(**plot_kwargs.get(target.name, {}))
    ax.grid(True)
    plt.gca().legend(handles=SOURCE_LEGEND_HANDLES, loc=(-0.15, 0.91))
    plot_last_update(y=-0.05)


def _plot_timeline_inset(start_date, end_date, x, y, dy=0.15, target_list=LTP):
    """
    """
    _start = date2num(numpy.datetime64(start_date))
    _end = date2num(numpy.datetime64(end_date))
    dx = (_end - _start) * 0.011
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
    plt.figure('IXPE LTP timeline', figsize=(18., 10.))
    LTP.plot()
    plt.gca().xaxis_date()
    plt.gca().yaxis.set_ticks([])
    plt.axis([None, None, -1., 1.])
    plt.tight_layout()
    _plot_timeline_inset('2022-01-31', '2022-02-25', 0.0075, 0.80)
    _plot_timeline_inset('2022-02-25', '2022-04-19', 0.0075, 0.05)
    _plot_timeline_inset('2022-05-11', '2022-07-08', 0.30, 0.62)
    _plot_timeline_inset('2022-08-27', '2022-10-03', 0.56, 0.25)
    plt.gca().legend(handles=SOURCE_LEGEND_HANDLES, loc=(0.20, 0.22))
    plot_last_update()



if __name__ == '__main__':
    print(TARGET_DICT)
    print(LTP)
    plot_mdp()
    plot_mdp(flux_range=False)
    plot_celestial()
    plot_galactic()
    plot_timeline()
    save_all_figures(IXPEOBSSIM_DATA)
    plt.show()

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


from enum import Enum, unique

from matplotlib.lines import Line2D

from ixpeobssim.utils.matplotlib_ import plt


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


# Collapse the enum into a list, which will be handy later.
TWG_LIST = list(TWG)

# A few matplotlib-related things, for plotting purposes.
_COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']
COLOR_DICT = {cls: _COLORS[i] for i, cls in enumerate(TWG_LIST)}
COLOR_DICT[None] = 'black'
SOURCE_LEGEND_HANDLES = []
for i, cls in enumerate(TWG):
    handle = Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_DICT[cls],
        markersize=10., label=cls.value)
    SOURCE_LEGEND_HANDLES.append(handle)



# This is compiled by hand and maintained to match the evolving ARTL.
#
# The dictionary is indexed by source name and the fields are (ra, dec, twg).
_TARGET_DATA = {
    '1ES 0229+200': (38.20256456, 20.2881926, TWG.BLAZAR),
    '1ES 1959+650': (299.99938562, 65.14851421, TWG.BLAZAR),
    '1RXS J170849.0': (257.20416667, -40.15277778, TWG.MAGNETAR),
    '3C 273': (187.27791535, 2.05238857, TWG.BLAZAR),
    '3C 279': (194.04652688, -5.7893144, TWG.BLAZAR),
    '3C 454.3': (343.49061658, 16.14821142, TWG.BLAZAR),
    '4U 0142+61': (26.59253, 61.75106, TWG.MAGNETAR),
    '4U 1626-67': (248.06996, -67.46091, TWG.AWDNS),
    '4U 1630-472': (248.506708, -47.393, TWG.ASMBH),
    '4U 1820-303': (275.919042, -30.361278, TWG.AWDNS),
    'BL Lac': (330.68038074, 42.27777178, TWG.BLAZAR),
    'Cas A': (350.85, 58.815, TWG.SNR),
    'Cen A': (201.36506288, -43.01911267, TWG.BLAZAR),
    'Cen X-3': (170.31288349, -60.62378542, TWG.AWDNS),
    'Circinus galaxy': (213.291458, -65.339222, TWG.RQAGN),
    'Crab': (83.63308333, 22.0145, TWG.PWN),
    'Cyg X-1': (299.59031591, 35.20160625, TWG.ASMBH),
    'Cyg X-2': (326.17147688, 38.32140718, TWG.AWDNS),
    'Cyg X-3': (308.107417, 40.95775, TWG.AWDNS),
    'GRO J1008-57': (152.4456843, -58.29322329, TWG.AWDNS),
    'GRS 1915+105': (288.798149, 10.945807, TWG.ASMBH),
    'GS 1826-238': (277.3675, -23.79694444, TWG.AWDNS),
    'GX 301-2': (186.65650295, -62.77034989, TWG.AWDNS),
    'GX 9+9': (262.93404167, -16.96136111, TWG.AWDNS),
    'Her X-1': (254.45754617, 35.34235762, TWG.AWDNS),
    'IC 4329A': (207.33025303, -30.3095067, TWG.RQAGN),
    'J0211+1051': (32.804905, 10.859666, TWG.BLAZAR),
    'LMC X-1': (84.91178513, -69.74320321, TWG.ASMBH),
    'MCG-5-23-16': (146.917319, -30.948734, TWG.RQAGN),
    'MSH 15-52' : (228.320833, -59.0816667, TWG.PWN),
    'Mrk 421': (166.113808, 38.20883287, TWG.BLAZAR),
    'Mrk 501': (253.46756952, 39.76016915, TWG.BLAZAR),
    'NGC 4151': (182.63573326, 39.40585098, TWG.RQAGN),
    'S5 0716+714': (110.47270204, 71.34343391, TWG.BLAZAR),
    'Sgr A complex': (266.41681662, -29.00782497, TWG.RQAGN),
    'SN 1006': (225.59208, -42.09694, TWG.SNR),
    'Tycho': (6.308299, 64.144305, TWG.SNR),
    'Vela Pulsar': (128.5, -45.83333333, TWG.PWN),
    'Vela X-1': (135.52858781, -40.55469345, TWG.AWDNS),
    'X Persei': (58.84615783, 31.04584604, TWG.AWDNS),
    'XTE J1701-462': (255.243583, -46.185722, TWG.AWDNS),
    'GRB 221009A': (288.254, 19.809, None),
    'EXO 2030+375': (308.06364577, 37.6374564, TWG.AWDNS),
    'PSR B0540-69': (85.04516, -69.33173, TWG.PWN),
}



# These were copied and formatted by hand from an excel spreadsheet maintained by
# Steve Odell. (At the cost of a small risk of editing errors, we preferred not
# read the excel file directly, as this would have added another dependence on
# ixoeobssim).
#
# The header row of the excel file reads:
#
# Source Name | (TWG() | Expos [d] | UID | Start | Stop | (Internal) | (Public) | Comments
#
# Note that we explicitly left out the fields for the TWG, as well as the
# internal and external data release.
_Y1_ARTL_DATA = (
    ('Cas A', 11.57, 1001301, '2022-01-11T11:45', '2022-01-29T12:34', ''),
    ('Cen X-3',	1.16, 1006501, '2022-01-29T13:10', '2022-01-31T07:10', ''),
    ('4U 0142+61',	6.94, 1003201, '2022-01-31T07:36', '2022-02-14T23:48', '1 of 2'),
    ('Cen A', 1.16, 1004301, '2022-02-15T00:38', '2022-02-17T13:33', ''),
    ('Her X-1', 2.31, 1001801, '2022-02-17T14:11', '2022-02-21T15:50', '1 of 2'),
    ('Crab', 0.58, 1001001, '2022-02-21T16:34', '2022-02-22T18:23', '1 of 2'),
    ('Her X-1', 1.16, 1001802, '2022-02-22T19:07', '2022-02-24T19:24', '2 of 2'),
    ('4U 0142+61',	1.97, 1003202,'2022-02-24T19:47', '2022-02-27T18:51', '2 of 2'),
    ('Sgr A complex', 4.05, 1003401, '2022-02-27T19:37', '2022-03-06T23:47', '1 of 2'),
    ('Crab', 0.58, 1001002, '2022-03-07T00:38', '2022-03-08T02:17', '2 of 2'),
    ('Mrk 501', 1.16, 1004501, '2022-03-08T02:59', '2022-03-10T08:30', ''),
    ('Sgr A complex', 7.52, 1003402, '2022-03-10T08:30', '2022-03-23T01:42', '2 of 2'),
    ('4U 1626-67', 2.31, 1002701, '2022-03-24T01:59', '2022-03-27T05:19', ''),
    ('Mrk 501', 1.16, 1004601, '2022-03-27T05:59', '2022-03-29T07:03', ''),
    ('GS 1826-238',	1.16, 1002801, '2022-03-29T07:25', '2022-03-31T08:57', ''),
    ('S5 0716+714',	4.63, 1005301, '2022-03-31T09:42', '2022-04-05T19:29', ''),
    ('Vela Pulsar', 5.79, 1001201, '2022-04-05T20:11', '2022-04-15T18:02', '1 of 2'),
    ('Vela X-1', 3.47, 1002501, '2022-04-15T18:11', '2022-04-21T12:17', ''),
    ('Vela Pulsar', 5.24, 1001202, '2022-04-21T12:24', '2022-04-30T10:02', '2 of 2'),
    ('Cyg X-2',	1.16, 1001601, '2022-04-30T11:04', '2022-05-02T11:07', ''),
    ('Cyg X-2', 0.57, 1006601, '2022-05-02T11:11', '2022-05-03T11:13', 'off-set'),
    ('1ES 1959+650', 0.58, 1006201, '2022-05-03T11:29', '2022-05-04T09:49', ''),
    ('Mrk 421',	1.16, 1003701, '2022-05-04T10:11', '2022-05-06T10:50', ''),
    ('BL Lac', 4.63, 1006301, '2022-05-06T11:30', '2022-05-14T12:27', ''),
    ('MCG-5-23-16', 0.58, 1003301, '2022-05-14T13:17', '2022-05-15T14:55', '1 of 2'),
    ('Cyg X-1',	3.47, 1002901, '2022-05-15T15:45', '2022-05-21T17:52', ''),
    ('MCG-5-23-16', 5.21, 1003302, '2022-05-21T18:42', '2022-05-31T03:47', '2 of 2'),
    ('3C 454.3', 1.16, 1005401, '2022-05-31T04:36', '2022-06-02T08:40', ''),
    ('3C 273', 1.16, 1005901, '2022-06-02T08:53', '2022-06-04T10:48', ''),
    ('Mrk 421',	1.16, 1003801, '2022-06-04T11:05', '2022-06-06T10:57', ''),
    ('1ES 1959+650', 0.58, 1006101, '2022-06-06T11:20', '2022-06-07T08:38', 'adjust'),
    ('Mrk 421',	1.16, 1003901, '2022-06-07T09:01', '2022-06-09T09:39', ''),
    ('1ES 1959+650', 2.46, 1006001, '2022-06-09T10:02', '2022-06-12T20:25', ''),
    ('3C 279', 3.24, 1005701, '2022-06-12T21:05', '2022-06-18T20:21', ''),
    ('Cyg X-1', 1.16, 1250101, '2022-06-18T21:03', '2022-06-21T21:03', 'ToO'),
    ('Tycho', 11.57, 1001401, '2022-06-21T21:22', '2022-07-07T00:00', ''),
    ('Cen X-3', 2.31, 1250201, '2022-07-04T06:22', '2022-07-07T12:44', 'ToO'),
    ('BL Lac', 1.44, 1006701, '2022-07-07T13:35', '2022-07-09T23:12', ''),
    ('Mrk 501', 1.16, 1004701, '2022-07-09T23:33', '2022-07-12T00:40', ''),
    ('Circinus galaxy', 9.26, 1003501, '2022-07-12T01:23', '2022-07-25T02:05', ''),
    ('SN 1006', 2.31, 1001501, '2022-07-25T02:20', '2022-07-29T02:19', 'Seg 1/2 Wrong coordinates'),
    ('GX 301-2', 3.47, 1002601, '2022-07-29T02:35', '2022-08-03T06:01', ''),
    ('SN 1006', 1.27, 1001502, '2022-08-03T06:16', '2022-08-05T14:18', 'Seg 2/2: Obs UID=1001599 Wrong coordinates'),
    ('SN 1006', 7.99, 1006801, '2022-08-05T14:21', '2022-08-19T12:32', ''),
    ('X Persei', 1.61, 1001701, '2022-08-19T13:22', '2022-08-22T05:39', 'Seg 1/2 interrupted by ToO'),
    ('4U 1630-472', 0.97, 1250301, '2022-08-22T06:28', '2022-08-23T22:37', 'Wrong coordinates'),
    ('4U 1630-472', 5.83, 1250401, '2022-08-23T22:47', '2022-09-02T18:54', 'ToO'),
    ('MSH 15-52', 10.48, 1001101, '2022-09-02T19:06', '2022-09-16T20:19', ''),
    ('X Persei', 1.28, 1001702, '2022-09-16T21:07', '2022-09-19T03:41', 'Seg 2/2 Obs UID=1001799 post ToO'),
    ('1RXS J170849.0', 0.08, 1003101, '2022-09-19T04:44', '2022-09-19T07:20', 'Seg 1/3'),
    ('1RXS J170849.0', 5.76, 1003102, '2022-09-19T17:30', '2022-09-29T12:01', 'Seg 2/3'),
    ('XTE J1701-462', 0.58, 1250601, '2022-09-29T12:09',	'2022-09-30T12:08', ''),
    ('1RXS J170849.0', 4.58, 1003103, '2022-09-30T12:14', '2022-10-08T11:19', 'Seg 3/3: Obs UID=1003199'),
    ('XTE J1701-462', 0.58, 1250602, '2022-10-08T11:27', '2022-10-09T11:25', ''),
    ('GX 9+9', 1.16, 1002401, '2022-10-09T11:41', '2022-10-11T13:13', ''),
    ('4U 1820-303', 0.12, 2002301, '2022-10-11T13:24', '2022-10-11T22:55', ''),
    ('GRB 221009A', 1.16, 2250101, '2022-10-11T23:15', '2022-10-14T00:46', 'ToO fast'),
    ('Cyg X-3',	3.23, 2001801, '2022-10-14T01:01', '2022-10-19T14:14',	'Seg 1/2'),
    ('LMC X-1',	6.94, 2001901, '2022-10-19T15:01', '2022-10-28T04:39', ''),
    ('1ES 1959+650', 2.31, 2004801, '2022-10-28T05:29', '2022-10-31T12:21', ''),
    ('Cyg X-3', 3.34, 2001802, '2022-10-31T12:35', '2022-11-06T09:34', 'Seg 2/2: Obs UID=2002899'),
    ('MCG-5-23-16', 3.71, 2003201, '2022-11-06T10:24', '2022-11-13T03:04', 'Seg 1/3'),
    ('GRO J1008-57', 0.98, 2003501, '2022-11-13T03:19', '2022-11-14T17:45', 'Re-released (gain issue)'),
    ('MCG-5-23-16', 1.91, 2003202, '2022-11-14T18:00', '2022-11-18T03:33', 'Seg 2/3'),
    ('GRO J1008-57', 1.23, 2003601, '2022-11-18T03:47', '2022-11-20T02:15', 'Late release (gain issue)'),
    ('MCG-5-23-16',	1.98, 2003203, '2022-11-20T02:30', '2022-11-23T18:29', 'Seg 3/3: Obs UID=2003299'),
    ('EXO 2030+375', 2.10, 2250201, '2022-11-23T19:19', '2022-11-27T13:10',	'ToO'),
    ('BL Lac', 1.89, 2005901, '2022-11-27T13:21', '2022-11-30T21:35', ''),
    ('Vela X-1', 3.19, 2005801, '2022-11-30T22:36', '2022-12-06T14:51', ''),
    ('Mrk 421',	0.87, 2004401, '2022-12-06T15:16', '2022-12-08T03:56', ''),
    ('NGC 4151', 7.63, 2003101, '2022-12-08T04:06', '2022-12-21T14:51', ''),
    ('Tycho', 2.57, 2001601, '2022-12-21T15:16', '2022-12-25T09:48', ''),
    ('Cyg X-3', 2.43, 2250301, '2022-12-25T10:05', '2022-12-29T17:45', 'ToO (slow)'),
    )

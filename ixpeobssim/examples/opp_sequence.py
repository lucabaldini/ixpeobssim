# Copyright (C) 2021, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public Licensese as published by
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

import os

from ixpeobssim import IXPEOBSSIM_DATA, IXPEOBSSIM_CALDB
import ixpeobssim.core.pipeline as pipeline
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.utils.time_ import days_to_seconds, string_to_met_utc


OBSERVING_PLAN = (
    ('crab_complex', '2022-02-19', days_to_seconds(1.)),
    ('dc1_cygx1', '2022-02-20', days_to_seconds(1.)),
    ('sgra', '2022-02-21', days_to_seconds(4.)),
    ('dc1_cygx1', '2022-02-25', days_to_seconds(1.)),
    ('opp_x_persei', '2022-02-26', days_to_seconds(3.)),
    ('opp_mrk_501', '2022-03-01', days_to_seconds(1.)),
    ('sgra', '2022-03-02', days_to_seconds(6.)),
    ('crab_complex', '2022-03-08', days_to_seconds(1.)),
)


DEFAULT_XPOBSSIM_KWARGS = dict(saa=True, occult=True, charging=True,
                               onorbitcalib=True, lv1a=True)

_base_path = os.path.join(IXPEOBSSIM_CALDB, 'bcf', 'chrgparams', 'ixpe_vanilla_d%d_chrgparams.fits')
CHRG_PARAMS_FILES = [_base_path % du_id for du_id in DU_IDS]


def test_times():
    """Test that the start time of each obervation coincides with the stop time
    of the previous one.
    """
    last_stop_met = None
    for obs_id, (model_name, start_date, duration) in enumerate(OBSERVING_PLAN):
        start_met = string_to_met_utc(start_date, lazy=True)
        stop_met = start_met + duration
        print('%s (%d): %.3f--%.3f' % (model_name, obs_id, start_met, stop_met))
        if last_stop_met is not None:
            assert start_met == last_stop_met
        last_stop_met = stop_met


def generate():
    """Generate a sequence of observations from one of the tentative versions
    of the observing plan for the first few weeks.
    """
    charge_maps = None
    for obs_id, (model_name, start_date, duration) in enumerate(OBSERVING_PLAN):
        outfile = os.path.join(IXPEOBSSIM_DATA, 'opp%02d_%s' % (obs_id, model_name))
        kwargs = dict(startdate=start_date, duration=duration, outfile=outfile)
        kwargs.update(chrgparams=CHRG_PARAMS_FILES)
        if charge_maps is not None:
            kwargs.update(chrgmaps=charge_maps)
        kwargs.update(**DEFAULT_XPOBSSIM_KWARGS)
        pipeline.setup(model=model_name)
        file_list = pipeline.xpobssim(**kwargs)
        charge_maps = pipeline.xpchrgmap(*file_list)


def run():
    """
    """
    generate()



if __name__ == '__main__':
    pipeline.bootstrap_pipeline(None)

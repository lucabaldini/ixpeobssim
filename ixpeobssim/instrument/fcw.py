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

"""On-orbit calibration facilities.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.instrument import DU_IDS, NUM_DETECTOR_UNITS
from ixpeobssim.utils.time_ import xTimeInterval


# pylint: disable=invalid-name


class xOnOrbitCalibrationRun(xTimeInterval):

    """Basic descriptor for a data-taking run with a calibration source from the FCW.
    """

    def __init__(self, start_met, stop_met, du_id, calibration_source):
        """Constructor.
        """
        xTimeInterval.__init__(self, start_met, stop_met)
        self.du_id = du_id
        self.calibration_source = calibration_source

    def __str__(self):
        """String formatting.
        """
        return 'Calibration run @ MET %s for DU %d' % (xTimeInterval.__str__(self), self.du_id)



class xOnOrbitCalibrationPattern(dict):

    """Convenience class representing the pattern for the on-orbit calibration.

    The basic idea, here, is that the detector units are calibrated one at a
    time, when the target is occoulted by the Earth and the observatory is not
    in the SAA. An optional demultiplier factor allows to calibrate each of
    the DUs every n orbits.

    The class is essentially a dictionary of xOnOrbitCalibrationRun objects,
    indexed by du_id---which is the typical way we retrieve them downstream.
    In addition, we keep a list of the run, preserving the time order, that
    is handy, e.g., when it comes to printing the object.
    """

    def __init__(self, octi_list, calibration_source, demultiplier=1):
        """Constructor.
        """
        dict.__init__(self, {du_id: [] for du_id in DU_IDS})
        self.calibration_source = calibration_source
        self.run_list = []
        for i, (start, stop) in enumerate(octi_list):
            if i % demultiplier == 0:
                du_id = (i // demultiplier) % NUM_DETECTOR_UNITS + 1
                run = xOnOrbitCalibrationRun(start, stop, du_id, calibration_source)
                self[du_id].append(run)
                self.run_list.append(run)

    def total_calibration_time(self, du_id):
        """Return the total calibration time for a given DU.
        """
        return sum([run.duration for run in self[du_id]])

    def num_calibration_runs(self, du_id):
        """Return the number of calibration runs for a given DU.
        """
        return len(self[du_id])

    def octi_data(self, du_id):
        """Return the onboard calibration time interval data to fill the proper
        extension in the output file.
        """
        runs = self[du_id]
        start = numpy.array([run.start_met for run in runs])
        stop = numpy.array([run.stop_met for run in runs])
        return (start, stop)

    def __str__(self):
        """String formatting.
        """
        text = 'List of calibration runs:'
        for i, run in enumerate(self.run_list):
            text += '\n[%3d] %s' % (i, run)
        for du_id in DU_IDS:
            text += '\nTotal calibration time for DU %d: %.3f s in %d run(s)' %\
                (du_id, self.total_calibration_time(du_id), self.num_calibration_runs(du_id))
        return text

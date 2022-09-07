#!/usr/bin/env python
#
# Copyright (C) 2020, the ixpeobssim team.
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

"""This is the simplest possible realization of the instrumental background,
where the energy spectrum is a power law fitted to the data from Bunner et al.
(`1978ApJ...220..261B <https://ui.adsabs.harvard.edu/abs/1978ApJ...220..261B/abstract>`_).
Here the authors provide the non X-ray background rates for their three detectors
and we are using values for the Neon-filled detector in Table 3 of the paper.

.. image:: figures/models/instrumental_bkg_spectrum.png

For completeness: this model component can be imported in any ixpeobssim
configuration file by simply adding the following line:

.. code-block:: python

    from ixpeobssim.config.instrumental_bkg import bkg

Then all you need to do is to add the bkg source to your region of interest.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.config import file_path_to_model_name
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.srcmodel.roi import xROIModel
from ixpeobssim.srcmodel.bkg import xPowerLawInstrumentalBkg


__model__ = file_path_to_model_name(__file__)

ra = 45.
dec = 45.


bkg = xPowerLawInstrumentalBkg()

ROI_MODEL = xROIModel(ra, dec, bkg)


def display(emin=1., emax=12.):
    """Display the source model.
    """
    energy = numpy.linspace(emin, emax, 100)
    # Energy spectrum
    plt.figure('%s spectrum' % __model__)
    plt.plot(energy, bkg.photon_spectrum(energy))
    setup_gca(xmin=emin, xmax=emax, ymin=bkg.photon_spectrum(emax),
              logx=True, logy=True, grids=True, **fmtaxis.spec)



if __name__ == '__main__':
    from ixpeobssim.config import bootstrap_display
    bootstrap_display()

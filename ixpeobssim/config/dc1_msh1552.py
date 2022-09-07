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

"""


Morphological and Spectral Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The basic morphological and spectral model for the nebula comes from a Chandra
event list, that we feed into ``ixpeobssim`` to be folded with the IXPE
response functions.

.. _figure-dc1_msh1552_chandra_image:
.. figure:: figures/dc1/dc1_msh1552_chandra_image.*
   :width: 80%

We cut out a small region at the position of the pulsar to avoid any problem with
pile up in Chandra and plug in our own pulsar model.


Polarization Model
~~~~~~~~~~~~~~~~~~

The baseline polarization model was created by Niccolo Bucciantini.

.. _figure-dc1_msh1552_polarization_degree:
.. figure:: figures/dc1/dc1_msh1552_polarization_degree.*
   :width: 80%

.. _figure-dc1_msh1552_polarization_angle:
.. figure:: figures/dc1/dc1_msh1552_polarization_angle.*
   :width: 80%

.. _figure-dc1_msh1552_polarization_model:
.. figure:: figures/dc1/dc1_msh1552_polarization_model.*
   :width: 80%




The PSR B1509-58 Pulsar
~~~~~~~~~~~~~~~~~~~~~~~

The pulsar is located at R.A. = 228.48133333, Dec. = -59.13577778. The basic
input for the modeling comes from the NuSTAR observation
`Ge Chen at al., 2015 <https://arxiv.org/abs/1507.08977v1>`_, where a full
phase-resolved spectral analysis is available.

We use a simple power-law parametrization, with a purely phenomenological
phase-dependence of the index, as illustrated in the figure.

.. _figure-dc1_msh1552_pulsar_index:
.. figure:: figures/dc1/dc1_msh1552_pulsar_index.*
   :width: 80%


The ephemeris are taken from
`Livingstone and Kaspi, 2011 <https://iopscience.iop.org/article/10.1088/0004-637X/742/1/31>`_:

.. code-block::

   Dates (Modified Julian Day) 	50148.096–55521.082
   Epoch (Modified Julian Day) 	52834.589000
   nu                           6.611515243850 s-1
   nudot0                       −6.694371307e-11 s-2
   nuddot                       1.9185594e-21 s-3


The polarization degree for the pulsar is taken as 35%, independently from the
pulse phase, with a polarization degree vs. pulse-phase profile provided by
Roger:

.. _figure-dc1_msh1552_pulsar_pol_ang:
.. figure:: figures/dc1/dc1_msh1552_pulsar_pol_ang.*
   :width: 80%


"""

import os

import numpy

from ixpeobssim import IXPEOBSSIM_CONFIG_FITS, IXPEOBSSIM_CONFIG_REG, IXPEOBSSIM_CONFIG_ASCII
from ixpeobssim.config import file_path_to_model_name, bootstrap_display
from ixpeobssim.config.dc1_bkg import instrumental_bkg
from ixpeobssim.config.dc1_utils import align_ephemeris
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.srcmodel.bkg import xGalacticBkg, xExtragalacticBkg
from ixpeobssim.srcmodel.img import xFITSImage
from ixpeobssim.srcmodel.polarization import constant, xStokesSkyMap
from ixpeobssim.srcmodel.roi import xChandraROIModel, xPeriodicPointSource, xChandraObservation
from ixpeobssim.srcmodel.spectrum import power_law
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.astro import read_ds9


#pylint: disable=invalid-name, unused-argument, redefined-outer-name


__model__ = file_path_to_model_name(__file__)


# Source coordinates, in decimal degrees.
SRC_NAME = 'MSH 15-52'
SRC_RA, SRC_DEC = 228.32083333, -59.08166667
SRC_L, SRC_B = 320.27811097, -1.07280689
PSR_NAME = 'PSR B1509-58'
PSR_RA, PSR_DEC = 228.48133333, -59.13577778
PSR_L, PSR_B = 320.32063283, -1.16174923

# Pointing coordinates
PNT_RA, PNT_DEC = PSR_RA, PSR_DEC

# Observation start date
START_DATE = '2022-02-01'
DURATION = 3000000.

# Interstellar absorption.
NH = 0.95e22

# Galactic diffuse.
GAL_BKG_RATE = 70.

# Phase-averaged pulsar integral flux.
# (Straight from table 3 of https://arxiv.org/pdf/1507.08977v1.pdf)
PSR_AVE_INTEGRAL_FLUX = 1.198
PSR_AVE_INT_FLUX_EMIN = 3.
PSR_AVE_INT_FLUX_EMAX = 79.
PSR_AVE_PL_NORM = 3.70e-3
PSR_AVE_PL_INDEX = 1.386

# Align the ephemeris to the start date.
PSR_EPHEMERIS = align_ephemeris(START_DATE, 52834.589000, 6.611515243850,
                                -6.694371307e-11, 1.9185594e-21)

# Constant polarization degree, as suggested by Roger.
PSR_POL_DEG = constant(0.35)

#Event file taken from the chandra catalog---the source name is listed as G320.4-1.2
# and not as MSH15-52:
# http://hea-www.cfa.harvard.edu/ChandraSNR/G320.4-01.2/ observation id is 754.
#X-ray size of the SNR is 18x11 arcminutes
EVENT_LIST_PATH = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'msh1552_acis_evt2_filtered.fits')

# This comes from the Chandra db?
IMG_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'msh1552_chandra_E2100-10000_FLUXED.fits')

# Created by Niccolo.
POL_DEG_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'MSH_Rad_PF.fits')
POL_ANG_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG_FITS, 'MSH_Rad_PA.fits')

# Region file for the PSR.
PSR_REG_FILE_PATH = os.path.join(IXPEOBSSIM_CONFIG_REG, 'msh1552_pulsar_position.reg')



def _pulse_profile_offset(r, pulsed_fraction):
    """Calculate the offset that it is necessary to add to a given pulse profile
    to achieve a target pulsed fraction.

    See the PF_area definition in
    https://iopscience.iop.org/article/10.3847/0004-637X/817/2/93
    """
    n = len(r)
    return (r.sum() * (1. / pulsed_fraction - 1.) - n * r.min() / pulsed_fraction) / n


def _parse_psr_pulse_profile(pulsed_fraction=0.7):
    """Parse the PSR B1509-58 spectral data from a csv file.

    The data were extracted from figure 1 of the paper, using the bottom pulse
    profile (3--79) keV, under the assumption that the pulse profile does not
    change with energy.

    Note we don't use the phase value in the file, as the underlying data points
    have 64 equi-spaced bins, and we can calculate the bin centers without
    relying on the manual image processing.
    """
    file_path = os.path.join(IXPEOBSSIM_CONFIG_ASCII, 'psr_b1509_58_pp.csv')
    _, y = numpy.loadtxt(file_path, delimiter=',', unpack=True)
    num_bins = 64
    x = numpy.linspace(0., 1. - 1. / num_bins, num_bins) + 0.5 / num_bins
    if pulsed_fraction is not None:
        y += _pulse_profile_offset(y, pulsed_fraction)
    return x, y


def _parse_psr_spectral_data():
    """Parse the PSR B1509-58 spectral data from a text file.

    The file (corresponding to Table 5 in the paper) was downloaded as ASCII from
    https://iopscience.iop.org/article/10.3847/0004-637X/817/2/93
    """

    def _parse_line(line):
        """Small nested function to parse the content of a single line in the
        input file.
        """
        data = line.split()
        min_phase, max_phase = [float(val) for val in data[0].split('-')]
        pl_index = float(data[1])
        pl_index_err = float(data[3])
        alpha = float(data[7])
        alpha_err = float(data[9])
        beta = float(data[10])
        beta_err = float(data[12])
        return numpy.array([min_phase, max_phase, pl_index, pl_index_err, alpha,
                            alpha_err, beta, beta_err])

    file_path = os.path.join(IXPEOBSSIM_CONFIG_ASCII, 'psr_b1509_58_spec_data.txt')
    data = None
    with open(file_path) as input_file:
        for line in input_file:
            while not line.startswith('0'):
                line = next(input_file)
            if data is None:
                data = _parse_line(line)
            else:
                data = numpy.vstack((data, _parse_line(line)))
    return data


def _parse_psr_pol_angle():
    """Parse the pulsar polarization angle as a function of the pulse phase.

    Data points courtesy of Roger Romani, see email on April 14, 2021.
    For completeness, he suggested to use a constant polarization degree of ~0.35.
    """
    file_path = os.path.join(IXPEOBSSIM_CONFIG_ASCII, 'psr_b1509_58_pol_data.txt')
    phase, pol_ang = numpy.loadtxt(file_path, unpack=True)
    # Add a point at phase = 1 to ensure continuity.
    phase = numpy.append(phase, 1.)
    pol_ang = numpy.append(pol_ang, pol_ang[0] - 180.)
    pol_ang = numpy.radians(pol_ang)
    return phase, pol_ang


# Parse the pulsar pulse profile---note at this point the pulse profile is in
# arbitrary units.
pp_phase, pp_norm = _parse_psr_pulse_profile()
# Scale to the phase-averaged PL normalization.
pp_norm = pp_norm / pp_norm.mean() * PSR_AVE_PL_NORM
psr_pp = xInterpolatedUnivariateSpline(pp_phase, pp_norm)

# Parse the pulsar spectral data.
min_phase, max_phase, pl_index, pl_index_err, alpha, alpha_err, beta, beta_err =\
    _parse_psr_spectral_data().T
phase = 0.5 * (min_phase + max_phase)
phase_err = 0.5 * (max_phase - min_phase)

# Parse the pulsar polarization data.
psr_pol_ang_spline = xInterpolatedUnivariateSpline(*_parse_psr_pol_angle(),
                                                   **fmtaxis.pp_pol_ang_rad)

# Pulsar phase-averaged spectrum
psr_ave_spec = power_law(PSR_AVE_PL_NORM, PSR_AVE_PL_INDEX)

def psr_pl_index(t):
    """Power-law index for the pulsar as a function of the pulse phase.

    This was done by hand, superimposing the line on the data point (i.e., it is
    not even a fit.)
    """
    return 1.34 +  1.5 * (abs(t - 0.38))**1.6

def psr_pl_norm(t):
    """Power-law normalization for the pulsar as a function of the pulse phase.
    """
    return psr_pp(t)

def psr_pol_ang(E, t, ra=None, dec=None):
    """Polarization angle for the pulsar as a function of the pulse phase.
    """
    return psr_pol_ang_spline(t)

# Pulsar phase-resolved spectrum.
psr_spec = power_law(psr_pl_norm, psr_pl_index)

# And, now, set up the nebula.
nebula_pol_map = xStokesSkyMap.load_from_pda(POL_DEG_FILE_PATH, POL_ANG_FILE_PATH)
nebula_pol_deg = nebula_pol_map.polarization_degree_model()
nebula_pol_ang = nebula_pol_map.polarization_angle_model()

# Define the actual ROI model.
ROI_MODEL = xChandraROIModel(EVENT_LIST_PATH, acis='I')
# Remove the pulsar from the source
psr_cut_reg = read_ds9(PSR_REG_FILE_PATH)[0]
psr_cookie_cut = xChandraObservation('pulsar', constant(0.), constant(0.),
                                     psr_cut_reg, exclude=True)
ROI_MODEL.add_source(psr_cookie_cut)
nebula = xChandraObservation(SRC_NAME, nebula_pol_deg, nebula_pol_ang)
ROI_MODEL.add_source(nebula)
psr = xPeriodicPointSource(PSR_NAME, PSR_RA, PSR_DEC, psr_spec, PSR_POL_DEG,
                           psr_pol_ang, PSR_EPHEMERIS, NH)
egb = xExtragalacticBkg(PNT_RA, PNT_DEC)
dge = xGalacticBkg(PNT_RA, PNT_DEC, GAL_BKG_RATE)
ROI_MODEL.add_sources(psr, egb, dge, instrumental_bkg)


def display_nebula():
    """Display the nebula.
    """
    plt.figure(f'{__model__} Chandra image')
    img = xFITSImage(IMG_FILE_PATH)
    img.plot(vmin=1e-7, stretch='log')
    plt.figure(f'{__model__} polarization degree')
    nebula_pol_map.plot_polarization_degree()
    plt.figure(f'{__model__} polarization angle')
    nebula_pol_map.plot_polarization_angle()
    plt.figure(f'{__model__} polarization model')
    img.plot(vmin=1e-7, stretch='log')
    nebula_pol_map.plot_arrows(nside=30, scale=20.)


def display_pulsar_model_building(log_parabola=False):
    """Display the ingredients going into building the pulsar model.
    """
    plt.figure(f'{__model__} pulsar pulse profile')
    plt.errorbar(pp_phase, pp_norm, fmt='o')
    setup_gca(xlabel='Pulse phase', ylabel='Pulse profile [a.u.]', xmin=0.,
              xmax=1., ymin=0., grids=True)
    plt.figure(f'{__model__} pulsar PL normalization')
    x = numpy.linspace(0., 1., 100)
    plt.plot(x, psr_pl_norm(x))
    setup_gca(xlabel='Pulse phase', ylabel='PL normalization', xmin=0., xmax=1., grids=True)
    plt.figure(f'{__model__} pulsar PL index')
    plt.errorbar(phase, pl_index, pl_index_err, phase_err, fmt='o')
    x = numpy.linspace(0., 1., 100)
    plt.plot(x, psr_pl_index(x))
    setup_gca(xlabel='Pulse phase', ylabel='PL index', xmin=0., xmax=1., grids=True)
    if log_parabola:
        plt.figure(f'{__model__} pulsar alpha')
        plt.errorbar(phase, alpha, alpha_err, phase_err, fmt='o')
        setup_gca(xlabel='Pulse phase', ylabel='$\\alpha$', xmin=0., xmax=1., grids=True)
        plt.figure(f'{__model__} pulsar beta')
        plt.errorbar(phase, beta, beta_err, phase_err, fmt='o')
        setup_gca(xlabel='Pulse phase', ylabel='$\\beta$', xmin=0., xmax=1., grids=True)
    plt.figure(f'{__model__} pulsar polarization angle')
    psr_pol_ang_spline.plot(overlay=True)
    setup_gca(grids=True)



def display():
    """Display the source model.
    """
    display_nebula()
    display_pulsar_model_building()



if __name__ == '__main__':
    print(psr.ephemeris_info())
    bootstrap_display()

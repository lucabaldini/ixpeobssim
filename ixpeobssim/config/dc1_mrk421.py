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
The input model for Markarian 421 comes with two independent components that we
shall refer to as `Faint` and `Bright`, respectively.

Spectral Model
~~~~~~~~~~~~~~

The two spectral components are modeled as simple power laws with  a normalization
that varies in time and fixed spectral indices.

The `Faint` component features a moderate secular variation of the order of 10%
of the average values, randomly generated with an interpolated cubic spline.

The `Bright` component has a negligible baseline, and a prominent peak
(shaped like a gamma function) more or less at the center of the observation.
Its spectral index is significantly harger than that of the `Faint` component.


.. _figure-dc1_mrk421_pl_normalization:
.. figure:: figures/dc1/dc1_mrk421_pl_normalization.*
   :width: 80%

.. _figure-dc1_mrk421_pl_index:
.. figure:: figures/dc1/dc1_mrk421_pl_index.*
   :width: 80%

   Power-law normalization and index for the two model components as a function
   of time, for the target observation period. (The gray vertical lines indicate
   the six observations.)

Below is the integral flux in the IXPE band for the two components and for their
sum.

.. _figure-dc1_mrk421_integral_flux:
.. figure:: figures/dc1/dc1_mrk421_integral_flux.*
   :width: 80%

   Integral flux for the two model components as a function of time, for the
   target observation period. (The gray vertical lines indicate
   the six observations.)


Polarization Model
~~~~~~~~~~~~~~~~~~

Both components feature a polarization degree increasing with energy
over the IXPE band---roughly speaking the polarization degree doubles
moving from 3 to 10 keV. The polarization vectors for the two components
are orthogonal to each other.

.. _figure-dc1_mrk421_pol_deg_energy:
.. figure:: figures/dc1/dc1_mrk421_pol_deg_energy.*
   :width: 80%

   Polarization degree as a function of the energy for the two coponents.


Although neither polarization component has an explicit time dependence (both
the polarization degree and polarization angle do not depend on time), the
fact that they are orthogonal generates interesting effects as the corresponding
spectra evolve with time. The specific profile of the `Bright` component and
the actual observation intervals have been admittedly cherry picked to
maximize the variety of combination we observe.

.. _figure-dc1_mrk421_pol_deg_time:
.. figure:: figures/dc1/dc1_mrk421_pol_deg_time.*
   :width: 80%

.. _figure-dc1_mrk421_pol_ang_time:
.. figure:: figures/dc1/dc1_mrk421_pol_ang_time.*
   :width: 80%

   Polarization degree and angle at a fixed energy (3 keV) as a function of
   time. (Note this is representative of the polarization at the peak of the
   sensitivity, but it is only a proxy of the broadband polarization we will
   measure in real life.)

"""


import numpy
import scipy.stats

from ixpeobssim.config import file_path_to_model_name, bootstrap_display
from ixpeobssim.config.dc1_bkg import instrumental_bkg
from ixpeobssim.config.dc1_utils import sum_pol_deg_ang

from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.srcmodel.bkg import xExtragalacticBkg, xGalacticBkg
from ixpeobssim.srcmodel.polarization import constant
from ixpeobssim.srcmodel.roi import xPointSource, xROIModel
from ixpeobssim.srcmodel.spectrum import power_law, pl_integral_flux
from ixpeobssim.utils.fmtaxis import fmtaxis
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, metplot
from ixpeobssim.utils.time_ import string_to_met_utc, met_to_num, met_to_mjd


#pylint: disable=invalid-name, unused-argument, too-many-locals


__model__ = file_path_to_model_name(__file__)


# Source coordinates, in decimal degrees.
SRC_NAME = 'Mrk 421'
SRC_RA, SRC_DEC = 166.113808, 38.20883287
SRC_L, SRC_B = 179.831675, 65.03148037

# Pointing coordinates
PNT_RA, PNT_DEC = SRC_RA, SRC_DEC

# Observation start dates and duration.
START_DATES = ('2022-04-23', '2022-05-07', '2022-05-21', '2022-05-28', '2022-06-04', '2022-06-19')
DURATION = 90000.
OBSERVABILITY_START_DATE = '2022-04-22T15:08:03.214428'
OBSERVABILITY_END_DATE = '2022-06-21T23:41:05.086307'
OBSERVABILITY_START_MET = string_to_met_utc(OBSERVABILITY_START_DATE, lazy=True)
OBSERVABILITY_END_MET = string_to_met_utc(OBSERVABILITY_END_DATE, lazy=True)

# Interstellar absorption.
NH = 1.45e20

# Galactic diffuse.
GAL_BKG_RATE = 10.

# Model parameters
FAINT_PL_NORM_MEAN = 0.03
FAINT_PL_NORM_RMS = 0.003
FAINT_PL_INDEX = 2.9
FAINT_POL_DEG_MIN = 0.15
FAINT_POL_DEG_MAX = 0.30
FAINT_POL_ANG = 42.

BRIGHT_PL_NORM_MIN = 0.0001
BRIGHT_PL_NORM_MAX = 0.1
BRIGHT_PEAK_TIME = string_to_met_utc('2022-05-15', lazy=True)
BRIGHT_PEAK_WIDTH = 250000.
BRIGHT_PEAK_ALPHA = 3.
BRIGHT_PL_INDEX = 2.5
BRIGHT_POL_DEG_MIN = 0.45
BRIGHT_POL_DEG_MAX = 0.90
BRIGHT_POL_ANG = FAINT_POL_ANG - 90.


def erf_pol_deg(min_, max_, center=6., width=2.):
    """Convenience function build the polarization degree as a function of time.
    """
    #pylint: disable=no-member
    def func(E, t=None, ra=None, dec=None):
        z = (E - center) / width
        return min_ + 0.5 * (max_ - min_) * (1. + scipy.special.erf(z / numpy.sqrt(2.)))
    return func


_met = numpy.linspace(OBSERVABILITY_START_MET, OBSERVABILITY_END_MET, 10)
_faint_pl_norm = FAINT_PL_NORM_MEAN + numpy.random.normal(0., FAINT_PL_NORM_RMS, size=_met.shape)
faint_pl_norm_spline = xInterpolatedUnivariateSpline(_met, _faint_pl_norm, k=3)

def faint_pl_norm(t):
    """Power-law normalization as a function of time for the faint component.
    """
    return faint_pl_norm_spline(t)

def faint_pl_index(t):
    """Power-law index as a function of time for the faint component.
    """
    if isinstance(t, numpy.ndarray):
        return numpy.full(t.shape, FAINT_PL_INDEX)
    return FAINT_PL_INDEX

faint_spec = power_law(faint_pl_norm, faint_pl_index)
faint_pol_deg = erf_pol_deg(FAINT_POL_DEG_MIN, FAINT_POL_DEG_MAX)
faint_pol_ang = constant(numpy.radians(FAINT_POL_ANG))


_met = numpy.linspace(OBSERVABILITY_START_MET, OBSERVABILITY_END_MET, 500)
peak = scipy.stats.gamma.pdf(_met, BRIGHT_PEAK_ALPHA, BRIGHT_PEAK_TIME, BRIGHT_PEAK_WIDTH)
peak /= peak.max()
_bright_pl_norm = BRIGHT_PL_NORM_MIN + (BRIGHT_PL_NORM_MAX - BRIGHT_PL_NORM_MIN) * peak
bright_pl_norm_spline = xInterpolatedUnivariateSpline(_met, _bright_pl_norm, k=3)

def bright_pl_norm(t):
    """Power-law normalization as a function of time for the bright component.
    """
    return bright_pl_norm_spline(t)

def bright_pl_index(t):
    """Power-law index as a function of time for the bright component.
    """
    if isinstance(t, numpy.ndarray):
        return numpy.full(t.shape, BRIGHT_PL_INDEX)
    return BRIGHT_PL_INDEX

bright_spec = power_law(bright_pl_norm, bright_pl_index)
bright_pol_deg = erf_pol_deg(BRIGHT_POL_DEG_MIN, BRIGHT_POL_DEG_MAX)
bright_pol_ang = constant(numpy.radians(BRIGHT_POL_ANG))

faint = xPointSource('Faint', SRC_RA, SRC_DEC, faint_spec, faint_pol_deg,
                     faint_pol_ang, column_density=NH)
bright = xPointSource('Bright', SRC_RA, SRC_DEC, bright_spec, bright_pol_deg,
                      bright_pol_ang, column_density=NH)

# Create the background components.
egb = xExtragalacticBkg(PNT_RA, PNT_DEC)
dge = xGalacticBkg(PNT_RA, PNT_DEC, GAL_BKG_RATE)


ROI_MODEL = xROIModel(PNT_RA, PNT_DEC, faint, bright, egb, dge, instrumental_bkg)


def display_observations(**kwargs):
    """Draw vertical lines and text labels for the source observations.
    """
    kwargs.setdefault('ls', 'dashed')
    kwargs.setdefault('lw', 1.)
    kwargs.setdefault('color', 'gray')
    for i, start_date in enumerate(START_DATES):
        obs_id = i + 1
        start = string_to_met_utc(start_date, lazy=True)
        end = start + DURATION
        plt.axvline(met_to_num(start), **kwargs)
        plt.axvline(met_to_num(end), **kwargs)
        ymin, ymax = plt.ylim()
        y0 = ymin + 0.75 * (ymax - ymin)
        plt.text(met_to_num(end) + 0.25, y0, f'Observation {obs_id}',
                 rotation=90., color=kwargs.get('color'))


def display_light_curve():
    """Display the light curve.
    """
    met = numpy.linspace(OBSERVABILITY_START_MET, OBSERVABILITY_END_MET, 500)
    faint_norm = faint_pl_norm(met)
    faint_index = faint_pl_index(met)
    bright_norm = bright_pl_norm(met)
    bright_index = bright_pl_index(met)
    emin, emax = 2., 8.
    faint_flux = pl_integral_flux(faint_norm, faint_index, emin, emax)
    bright_flux = pl_integral_flux(bright_norm, bright_index, emin, emax)

    plt.figure(f'{__model__} pl normalization')
    metplot(met, faint_norm, label='Faint')
    metplot(met, bright_norm, label='Bright')
    metplot(met, faint_norm + bright_norm, label='Total')
    setup_gca(ylabel='Power-law normalization at 1 keV [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]',
              xlabel='Date', grids=True, legend=True)
    display_observations()

    plt.figure(f'{__model__} pl index')
    metplot(met, faint_index, label='Faint')
    metplot(met, bright_index, label='Bright')
    setup_gca(ylabel='Power-law index', xlabel='Date', grids=True, legend=True,
              ymin=2., ymax=3.)
    display_observations()

    plt.figure(f'{__model__} integral flux')
    metplot(met, faint_flux, label='Faint')
    metplot(met, bright_flux, label='Bright')
    metplot(met, faint_flux + bright_flux, label='Total')
    setup_gca(ylabel='Integral flux %d--%d keV [erg cm$^{-2}$ s$^{-1}$]' % (emin, emax),
              xlabel='Date', grids=True, legend=True)
    display_observations()


def display_polarization():
    """Display the polarization.
    """
    energy = numpy.linspace(1., 12., 100)
    plt.figure(f'{__model__} pol deg energy')
    plt.plot(energy, faint_pol_deg(energy), label='Faint')
    plt.plot(energy, bright_pol_deg(energy), label='Bright')
    setup_gca(xmin=energy.min(), xmax=energy.max(), ymin=0., grids=True,
              legend=True, **fmtaxis.ene_pol_deg)

    met = numpy.linspace(OBSERVABILITY_START_MET, OBSERVABILITY_END_MET, 500)
    ref_energy = 3.
    E = numpy.full(met.shape, ref_energy)
    faint_flux = faint_pl_norm(met) * E**faint_pl_index(met)
    faint_pd = faint_pol_deg(E, met)
    faint_pa = faint_pol_ang(E, met)
    bright_flux = bright_pl_norm(met) * E**bright_pl_index(met)
    bright_pd = bright_pol_deg(E, met)
    bright_pa = bright_pol_ang(E, met)
    _, pol_deg, pol_ang = sum_pol_deg_ang((faint_flux, faint_pd, faint_pa),
                                          (bright_flux, bright_pd, bright_pa))
    plt.figure(f'{__model__} pol deg time')
    metplot(met, faint_pd, label='Faint')
    metplot(met, bright_pd, label='Bright')
    metplot(met, pol_deg, label='Total')
    setup_gca(ylabel='Polarization degree @%.2f keV' % ref_energy, xlabel='Date',
              ymin=0., ymax=0.5, grids=True, legend=True)
    display_observations()

    plt.figure(f'{__model__} pol ang time')
    metplot(met, numpy.degrees(faint_pa), label='Faint')
    metplot(met, numpy.degrees(bright_pa), label='Bright')
    metplot(met, numpy.degrees(pol_ang), label='Total')
    setup_gca(ylabel='Polarization angle @%.2f keV [deg]' % ref_energy,
              xlabel='Date', grids=True, legend=True)
    display_observations()


def display_optical_polarization():
    """Display the optical polarization.
    """
    pol_deg_err = 0.20
    pol_ang_err = 3.0
    rmag_err = 0.05
    mjd = []
    pol_ang = []
    pol_deg = []
    rmag = []
    for i, start_date in enumerate(START_DATES):
        met = string_to_met_utc(start_date, lazy=True)
        _mjd = met_to_mjd(met)
        _mjd += numpy.random.uniform(-2, 2.)
        mjd.append(_mjd)
        if i == 2:
            pd = 6. + numpy.random.normal(0, pol_deg_err)
            pa = BRIGHT_POL_ANG + numpy.random.normal(0, pol_ang_err)
            r = 12.0 + numpy.random.normal(0, rmag_err)
        else:
            pd = 3. + numpy.random.normal(0, pol_deg_err)
            pa = FAINT_POL_ANG + numpy.random.normal(0, pol_ang_err)
            r = 12.5 + numpy.random.normal(0, rmag_err)
        print('MJD %.3f, pd = %.2f, pa = %.2f, R = %.2f' % (_mjd, pd, pa, r))
        pol_deg.append(pd)
        pol_ang.append(pa)
        rmag.append(r)
    pol_deg = numpy.array(pol_deg)
    pol_ang = numpy.array(pol_ang)
    rmag = numpy.array(rmag)
    fig, axs = plt.subplots(3, num=f'{__model__} optical polarization', sharex=True)
    axs[0].errorbar(mjd, rmag, rmag_err, fmt='o')
    plt.sca(axs[0])
    setup_gca(ylabel='R [mag]', grids=True, ymin=11.75, ymax=12.75)
    axs[1].errorbar(mjd, pol_deg, pol_deg_err, fmt='o')
    plt.sca(axs[1])
    setup_gca(ylabel='Pol. degree [%]', grids=True, ymin=0., ymax=9.)
    axs[2].errorbar(mjd, pol_ang, pol_ang_err, fmt='o')
    plt.sca(axs[2])
    setup_gca(xlabel='Time [MJD]', ylabel='Pol. angle [deg]', grids=True, ymin=-59., ymax=59.)
    fig.align_ylabels()


def display():
    """Main display() function.
    """
    display_light_curve()
    display_polarization()
    display_optical_polarization()



if __name__ == '__main__':
    bootstrap_display()

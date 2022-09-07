#!/urs/bin/env python
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


import os
import numpy

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from ixpeobssim import IXPEOBSSIM_SRCMODEL
from ixpeobssim.utils.logging_ import abort
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear



def mapped_column_density_HI(ra, dec, map_name='LAB'):
    """Return the mapped H_I column density at a given position in the sky.

    The value is read from one of the available input maps. Note that the
    data in the maps are stored in Galactic coordinates, while we're
    giving Celestial coordinates in input, here---the transformation is
    handled internally.

    Arguments
    ---------
    ra : float
        Right ascension of the source (in decimal degrees).

    dec: float
        Declination of the source (in decimal degrees).

    map: str
        The HI column density map to use. Can be either 'LAB' (LAB survey)
        or 'DL' (Dickey & Lockman).
    """
    # Make sure the map name makes sense.
    assert map_name in ['LAB', 'DL']
    # Transform from Celestial to Galactic coordinates.
    gal_coords = SkyCoord(ra, dec, unit='deg').galactic
    l, b = gal_coords.l.degree, gal_coords.b.degree
    # Open the selected input FITS map and grab the values.
    file_path = os.path.join(IXPEOBSSIM_SRCMODEL,'fits','h1_nh_%s.fits' %\
                             map_name)
    if not os.path.exists(file_path):
        abort('Could not find %s' % file_path)
    with fits.open(file_path) as hdu_list:
        _wcs = WCS(hdu_list[0].header)
        _data = hdu_list[0].data
    row, col = [int(item) for item in _wcs.wcs_world2pix(l, b, 1)]
    return _data[col, row]



class xInterstellarAbsorptionModel:

    """Class implementing the insterstellar absorption model using the
    Wisconsin (Morrison and McCammon; ApJ 270, 119) cross-sections.

    Here we essentially read the numbers in table 2 from the paper and
    build an interpolated univariate linear spline with the photoabsorption
    cross section values. The cross section is parametrized as a set of
    piecewise quadratic functions fitted to the data in energy bins---see the
    figure below.

    .. image:: ../docs/figures/gabs_xsection.png

    Example
    -------
    >>> from ixpeobssim.srcmodel.gabs import xInterstellarAbsorptionModel
    >>> model = xInterstellarAbsorptionModel()
    # Build a spline with the interstellar transmission factor as a function
    # of the photon energy
    >>> trans = model.transmission_factor(1.e22)
    >>> trans.plot()
    """

    def __init__(self, num_samples=1000):
        """Constructor.

        Arguments
        ---------
        num_samples : int
            The number of data points used to sample (logarithmically) the
            photon energies when tabulating the absorption cross section
            internally.
        """
        file_path = os.path.join(IXPEOBSSIM_SRCMODEL, 'ascii', 'XsecFitParam.txt')
        if not os.path.exists(file_path):
            abort('Could not find %s' % file_path)
        # Read the data from file and calculate the minimum and maximum
        # energies.
        e_lo, e_hi, c0, c1, c2 = numpy.loadtxt(file_path, delimiter=',',
                                               unpack=True)
        emin = e_lo[0]
        emax = e_hi[-1]
        # Sample the energy logarithmically between emin and emax.
        _x = numpy.logspace(numpy.log10(emin), numpy.log10(emax), num_samples)
        # Here comes the interesting part---the cross section is tabulated
        # by means of a set of piecewise quadratic functions, and the
        # three coefficients are provided in each energy bin. We do some
        # magic with the numpy.digitize() function, returning for each value
        # in _x its bin index in the original table. Note that we need to
        # clip the array removing negative numbers, as the double log/exp
        # operations done in defining the binning where driving the first
        # index to -1.
        _bin = (numpy.digitize(_x, e_lo) - 1).clip(min=0)
        # Calculate the cross section values---compact, isn't it?
        _y = 1.0e-24*(c0[_bin] + c1[_bin]*_x + c2[_bin]*(_x**2.))/(_x**3.)
        # And, finally, build the basic spline underlying the model.
        _fmt = dict(xlabel='Energy [keV]', ylabel='$\\sigma_{abs}$ [cm$^2$]')
        self.xsection = xInterpolatedUnivariateSplineLinear(_x, _y, **_fmt)

    def __call__(self, energy):
        """Return the value(s) of the photoelectric absorption cross section
        of the model for a value (or an array of values) of energy.

        Arguments
        ---------
        energy : float or array
            The energy (or energies) at which the cross section should be
            calculated.
        """
        return self.xsection(energy)

    def xsection_ecube(self):
        """Return a spline with the values of absorption cross section
        multiplied by the cube of the energy, as a function of the energy.

        This is essentially the same plot in Figure 1 of the original paper, as
        illustrated in the following figure

        .. image:: ../docs/figures/gabs_xsection_ecube.png
        """
        _x = self.xsection.x
        _y = self(_x)*(_x**3.)
        _fmt = dict(xlabel='Energy [keV]',
                    ylabel='$\\sigma_{abs} \\times E^3$ [cm$^2$~keV$^3$]')
        return xInterpolatedUnivariateSplineLinear(_x, _y, **_fmt)

    def transmission_factor(self, column_density):
        """Return the transmission factor for a given column density.

        This is essentially returning

        .. math::
            \\varepsilon = \\exp(-n_H\\sigma)

        .. image:: ../docs/figures/gabs_trans_samples.png

        Arguments
        ---------
        column_density : float
            The column density at which the transmission factor should be
            calculated.

        Warning
        -------
        We do have an issue with the extrapolation, here, as at the current
        stage there is no guarantee that the result of the spline evaluation
        would be <= 1. We could set the ext class member of the spline to 3
        before reurning it, but event that would not be right in general.
        This is probably not a huge issue, but it should be addressed
        properly.
        """
        _x = self.xsection.x
        _y = numpy.exp(-column_density*self(_x))
        _fmt = dict(xlabel='Energy [keV]', ylabel='Transmission factor')
        return xInterpolatedUnivariateSplineLinear(_x, _y, **_fmt)

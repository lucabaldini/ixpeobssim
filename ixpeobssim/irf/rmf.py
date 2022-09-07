# Copyright (C) 2015--2022, the ixpeobssim team.
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

"""Definition of the response matrix.
"""

from __future__ import print_function, division

import numpy

from ixpeobssim.core.fitting import fit
from ixpeobssim.core.modeling import xGaussian
from ixpeobssim.core.rand import xUnivariateAuxGenerator
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.irf.base import xResponseBase
from ixpeobssim.utils.matplotlib_ import plt, last_line_color, setup_gca

# pylint: disable=invalid-name, too-many-ancestors




class xEnergyDispersionMatrix(xUnivariateAuxGenerator):

    """Class encapsulating the energy dispersion matrix, as stored in the
    MATRIX extension of a .rmf file.

    We don't downsample the matrix to avoid the failure of the fit procedure in
    XSPEC.

    Arguments
    ---------
    hdu : FITS hdu
       The MATRIX hdu in the .rmf FITS file.
    """

    def __init__(self, hdu):
        """Constructor.
        """
        data = hdu.data
        _x = numpy.arange(0, len(data['MATRIX'][0]), 1)
        _y = 0.5 * (data['ENERG_LO'] + data['ENERG_HI'])
        _z = data['MATRIX'].T
        fmt = dict(xlabel='Channel', ylabel='Energy [keV]')
        xUnivariateAuxGenerator.__init__(self, _x, _y, _z, **fmt)

    @staticmethod
    def digitize(value, dtype=numpy.int16):
        """Digitize an energy value.

        This is essentially rounding an array to the nearest integer and
        cast the result to the specified type.
        """
        return numpy.ndarray.astype(numpy.rint(value), dtype)

    def rvs(self, aux):
        """Overloaded method.

        This is the original rvs() overloade method that was initially
        implemented in the energy dispersion code, and it is digitizing
        the measured energy before returning it, just like the detector
        would do in real life.

        We discovered along the way that it is handy to keep track of the
        measured (smneared) energy *before* the digitization (i.e., with all
        the significant digits), as this allows to apply corrections (e.g., for
        the GEM charging) without running into rounding problems. This is why
        we added the rvs2() class method, which is returning both.
        """
        val = xUnivariateAuxGenerator.rvs(self, aux)
        return self.digitize(val)

    def rvs2(self, aux):
        """Alterntative rvs() implementation where we return the measured
        energy both *before* and *after* the digitization.

        Note that the energy before the digitization is still in the pha space,
        and needs to be explicitely converted to physical units (keV), if
        that it what one is interested. When that's the case, the operation is
        deferred to the xEnergyDispersion class, who is aware of both the
        energy dispersion matric and the energy bounds.
        """
        val = xUnivariateAuxGenerator.rvs(self, aux)
        return val, self.digitize(val)



class xEnergyDispersionBounds(xInterpolatedUnivariateSplineLinear):

    """Class encapsulating the bounds for the energy dispersion matrix, as
    stored in the EBOUNDS extension of a .rmf file.

    This is essentially a univariate spline with the channel number on the
    x axis and the average energy of the corresponding interval on the y
    axis. At the same time we do store the bounds of the intervals themselves.
    """

    def __init__(self, hdu):
        """Constructor.
        """
        _bounds = hdu.data
        self.chans = _bounds['CHANNEL']
        self.emin = _bounds['E_MIN']
        self.emax = _bounds['E_MAX']
        _x = self.chans
        _y = 0.5 * (self.emin + self.emax)
        fmt = dict(xlabel='Channel', ylabel='Average energy [keV]')
        xInterpolatedUnivariateSplineLinear.__init__(self, _x, _y, **fmt)

    def num_channels(self):
        """Return the number of channels.
        """
        return len(self.chans)

    def min(self):
        """Return the minimum energy of the first channel.
        """
        return self.emin[0]

    def max(self):
        """Return the maximum energy of the last channel.
        """
        return self.emax[-1]

    def channel_to_energy(self, channel):
        """Convert a channel number to energy.
        """
        return self(channel)

    def energy_to_channel(self, energy):
        """Convert a physical energy (in keV) into a channel number.

        This is using a simple binary search on the underlying energy bounds.
        """
        if not isinstance(energy, numpy.ndarray):
            energy = numpy.array(energy)
        ch = numpy.searchsorted(self.emax, energy, side='left')
        if numpy.isscalar(ch):
            return numpy.int16(ch)
        return numpy.ndarray.astype(ch, numpy.int16)



class xEnergyDispersion(xResponseBase):

    """Class representing the energy dispersion.

    Arguments
    ---------
    file_path : str
        The path to the .rmf FITS file containing the energy dispersion tables.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xResponseBase.__init__(self, file_path, 'rmf')
        self.matrix = xEnergyDispersionMatrix(self.hdu_list['MATRIX'])
        self.ebounds = xEnergyDispersionBounds(self.hdu_list['EBOUNDS'])

    def channel_to_energy(self, channel):
        """Convenient proxy to the underlying xEnergyDispersionBounds method.
        """
        return self.ebounds.channel_to_energy(channel)

    def energy_to_channel(self, energy):
        """Convenient proxy to the underlying xEnergyDispersionBounds method.
        """
        return self.ebounds.energy_to_channel(energy)

    @staticmethod
    def pha_to_pi(pha):
        """Simple placeholder method to convert a PHA value to pulse invariant.

        Note
        ----
        Since we're not handling the gain corrections, yet, this is simply
        converting the values from integer to floating point numbers.
        """
        return pha.astype(numpy.double)

    def pha_analysis(self, energy):
        """Perform a pulse-height analysis on a given array of energies.
        """
        pha = self.energy_to_channel(energy)
        pi = self.pha_to_pi(pha)
        return pha, pi

    @staticmethod
    def scale_energy(energy):
        """Hook to allow to apply distorsions to the energy response of the
        detector.

        In this default implementation this is just returning the energy,
        untouched, but the method can be overridden to simulate a variety
        of systematic effects.
        """
        return energy

    def convolve_energy(self, mc_energy):
        """Convolve an array of mc_energy with the energy dispersion.

        This function is less trivial than the name might suggest, and it is
        performing several different operations at the same time, namely:

        * call xEnergyDispersionMatrix.rvs2() to smear the true energy and
          get the measured energy (in channel space) both before and after
          the digitization;
        * convert the measured energy before digitization to keV;
        * convert the pha to pulse invariant.
        """
        # Smear the Monte Carlo energy with the energy dispersion, and retrieve
        # the measured energy (in channel space) both before and after the
        # digitization.
        energy, pha = self.matrix.rvs2(self.scale_energy(mc_energy))
        # Convert the measured energy from channel to keV.
        energy = self.ebounds(energy)
        # Calculate the naive pulse invariant.
        pi = self.pha_to_pi(pha)
        return energy, pha, pi

    def plot(self):
        """Plot the energy dispersion.
        """
        plt.figure('%s energy dispersion' % self.base_name)
        self.matrix.plot()
        plt.figure('%s energy resolution' % self.base_name)
        for energy in [2.7, 5.9, 10.]:
            # Retrieve the energy dispersion 1-d spline at a given energy
            pdf = self.matrix.slice(energy)
            # Initial fit with a Gaussian...
            ch = self.energy_to_channel(energy)
            model = fit(xGaussian(), pdf.x, pdf.y, p0=(1., ch, 20.))
            # ... refine the fit around the peak.
            xmin = model.peak - 1.5 * model.sigma
            xmax = model.peak + 2.5 * model.sigma
            model = fit(model, pdf.x, pdf.y, xmin=xmin, xmax=xmax)
            # Retrieve the energy resolution.
            eres = 100. * model.resolution()
            pdf.plot()
            label = '$\\sigma_E$/E = %.1f%% FWHM @ %.2f keV' % (eres, energy)
            x = model.peak
            y = pdf(x)
            plt.text(x - 28, y, label, va='bottom', color=last_line_color())
        _, ymax = plt.ylim()
        setup_gca(ymin=0., ymax=1.15 * ymax)

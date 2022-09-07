#!/urs/bin/env python
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


from __future__ import print_function, division

import unittest
import sys
import os

import numpy
from astropy.io import fits

from ixpeobssim.utils.logging_ import logger
from ixpeobssim import IXPEOBSSIM_TEST_DATA, IXPEOBSSIM_CONFIG_ASCII, IXPEOBSSIM_AUXFILES
from ixpeobssim.srcmodel.magnetar import xMagnetarModelsT2020, parse_data_file,\
    xMagnetarTableModelT2020, xMagnetarTableModelT2020QedOff, parse_legacy_data
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
if sys.flags.interactive:
    plt.ion()


class TestMagnetarTableModels(unittest.TestCase):

    """Unit test for the magnetar table models.
    """

    @classmethod
    def setUpClass(cls):
        """
        """
        if not xMagnetarTableModelT2020.auxfiles_missing():
            cls.model = xMagnetarTableModelT2020()
        if not xMagnetarTableModelT2020QedOff.auxfiles_missing():
            cls.model_qedoff = xMagnetarTableModelT2020QedOff()

    def test_old_1708(self, integral_flux=2.43e-11):
        """Comparison with the old 1708 model. Note that the two don't need to be
        equal, as Roberto pointed out in one of his emails, as the underlying
        assumptions are different. Yet, this is still useful to check the
        orders of magnitude.
        """
        if xMagnetarTableModelT2020.auxfiles_missing():
            return
        file_name = 'axp_1rxs_j1708_phase_energy_1708-90-60-05-034.dat'
        file_path = os.path.join(IXPEOBSSIM_CONFIG_ASCII, file_name)
        phase, energy, flux, pol_deg, pol_ang = parse_legacy_data(file_path, False, 0.53, 9.5)
        params = (xMagnetarModelsT2020.BLACKBODY, 89.9, 60., 0.5, 0.34)
        spec_spline, pol_deg_spline, pol_ang_spline = self.model.interpolated_splines(*params, integral_flux)
        plt.figure('1708 Spectrum comparison')
        plt.plot(energy, flux[1], 'o', label='ASCII file')
        spec_spline.hslice(phase[1]).plot(label='Table model')
        setup_gca(xlabel='Energy [keV]', ylabel='Flux', grids=True, legend=True, logy=True)
        plt.figure('1708 Polarization degree comparison')
        plt.plot(energy, pol_deg[1], 'o', label='ASCII file')
        pol_deg_spline.hslice(phase[1]).plot(label='Table model')
        setup_gca(xlabel='Energy [keV]', ylabel='Polarization degree', grids=True, legend=True)
        plt.figure('1708 Polarization angle comparison')
        plt.plot(energy, pol_ang[1], 'o', label='ASCII file')
        pol_ang_spline.hslice(phase[1]).scale(180. / numpy.pi).plot(label='Table model')
        setup_gca(xlabel='Energy [keV]', ylabel='Polarization angle [$^\\circ$]', grids=True, legend=True)

    def test_base(self):
        """
        """
        if xMagnetarTableModelT2020.auxfiles_missing():
            return
        params = (xMagnetarModelsT2020.BLACKBODY, 89.9, 60., 0.5, 0.34)
        spec_spline, pol_deg_spline, pol_ang_spline = self.model.nearest_splines(*params, 3.e-11)
        plt.figure('Spectrum 2d')
        spec_spline.plot()
        plt.figure('Polarization degree 2d')
        pol_deg_spline.plot()
        plt.figure('Polarization angle 2d')
        pol_ang_spline.plot()

    def test_nearest(self):
        """
        """
        if xMagnetarTableModelT2020.auxfiles_missing():
            return
        file_path = os.path.join(IXPEOBSSIM_TEST_DATA, 'poldatatrans_89.90_60.00_05_034.dat')
        phase, energy, flux, pol_deg, pol_ang = parse_legacy_data(file_path, False, 0.53, 9.5)
        self.assertTrue(numpy.allclose(energy, self.model.energy))
        self.assertTrue(numpy.allclose(phase, self.model.phase))
        params = (xMagnetarModelsT2020.BLACKBODY, 89.9, 60., 0.5, 0.34)
        kwargs = dict(pad_phase=False, pad_energy=False)
        spec_spline, pol_deg_spline, pol_ang_spline = self.model.nearest_splines(*params, None, **kwargs)
        self.assertTrue(numpy.allclose(pol_deg, pol_deg_spline.z.T))
        self.assertTrue(numpy.allclose(pol_ang, numpy.degrees(pol_ang_spline.z.T)))
        plt.figure('Spectrum comparison')
        plt.plot(energy, flux[0], 'o', label='ASCII file')
        spec_spline.hslice(phase[0]).plot(label='Table model')
        setup_gca(xlabel='Energy [keV]', ylabel='Flux', grids=True, legend=True, logy=True)
        plt.figure('Polarization degree comparison')
        plt.plot(energy, pol_deg[0], 'o', label='ASCII file')
        pol_deg_spline.hslice(phase[0]).plot(label='Table model')
        setup_gca(xlabel='Energy [keV]', ylabel='Polarization degree', grids=True, legend=True)
        plt.figure('Polarization angle comparison')
        plt.plot(energy, pol_ang[0], 'o', label='ASCII file')
        pol_ang_spline.hslice(phase[0]).scale(180. / numpy.pi).plot(label='Table model')
        setup_gca(xlabel='Energy [keV]', ylabel='Polarization angle [$^\\circ$]', grids=True, legend=True)

    def test_interpolate_indices(self):
        """
        """
        if xMagnetarTableModelT2020.auxfiles_missing():
            return
        params = xMagnetarModelsT2020.BLACKBODY, 89.9, 60., 0.5, 0.34
        indices, weights = self.model._I.interpolate_indices(*params)
        self.assertTrue((indices == [1, 6, 4, 2, 2]).all())
        self.assertTrue(numpy.allclose(weights, [1.]))
        params = xMagnetarModelsT2020.BLACKBODY, 89.9, 60., 0.45, 0.34
        indices, weights = self.model._I.interpolate_indices(*params)
        self.assertTrue((indices == [[1, 6, 4, 1, 2], [1, 6, 4, 2, 2]]).all())
        self.assertTrue(numpy.allclose(weights, [0.5, 0.5]))

    def test_interpolate(self):
        """
        """
        if xMagnetarTableModelT2020.auxfiles_missing():
            return
        params = xMagnetarModelsT2020.BLACKBODY, 89.9, 60., 0.5, 0.34
        interp = self.model._interpolate_data(*params)
        nearest = self.model._nearest_data(*params)
        self.assertTrue(numpy.allclose(interp, nearest))
        for delta_phi in (0.5, 0.55, 0.6):
            params = xMagnetarModelsT2020.BLACKBODY, 89.9, 60., delta_phi, 0.34, None
            spec_spline, pol_deg_spline, pol_ang_spline = self.model.interpolated_splines(*params)
            plt.figure('Spectrum interpolation')
            spec_spline.hslice(0.5).plot(label='delta_phi = %.3f' % delta_phi)
            plt.legend()
            plt.figure('Polarization degree interpolation')
            pol_deg_spline.hslice(0.5).plot(label='delta_phi = %.3f' % delta_phi)
            plt.legend()
            plt.figure('Polarization angle interpolation')
            pol_ang_spline.hslice(0.5).plot(label='delta_phi = %.3f' % delta_phi)
            plt.legend()

    def test_continuity(self):
        """Test the the phase profiles are continuous (i.e., the answer at 0 is
        the same at that for phase = 1).
        """
        if xMagnetarTableModelT2020.auxfiles_missing():
            return
        params = (xMagnetarModelsT2020.BLACKBODY, 89.9, 60., 0.5, 0.34)
        integral_flux = 3.5e-11
        spec_spline, pol_deg_spline, pol_ang_spline = self.model.interpolated_splines(*params, integral_flux, 2., 10.)
        for energy in numpy.linspace(2., 8., 10):
            for spline in spec_spline, pol_deg_spline, pol_ang_spline:
                self.assertAlmostEqual(spline(energy, 0.), spline(energy, 1.))

    def test_normalization(self):
        """
        """
        if xMagnetarTableModelT2020.auxfiles_missing():
            return
        params = (xMagnetarModelsT2020.BLACKBODY, 89.9, 60., 0.5, 0.34)
        integral_flux = 3.5e-11
        spec_spline, _, _ = self.model.interpolated_splines(*params, integral_flux, 2., 10.)
        logger.info('Energy grid: %s', spec_spline.x)
        logger.info('Phase grid: %s', spec_spline.y)
        plt.figure('Test normalization')
        for phase in (0., 0.25, 0.5, 0.75, 1.):
            spec_spline.hslice(phase).plot(label='Phase = %.2f' % phase)
        plt.legend()

    def test_pol_deg_negative_values(self):
        """Test possible negative (i.e., unphysical) values of the polarization
        degree.
        """
        if xMagnetarTableModelT2020.auxfiles_missing():
            return
        params = (xMagnetarModelsT2020.SOLID_SURFACE_FIXED_IONS, 89.9, 60., 0.5, 0.34)
        integral_flux = 3.5e-11
        spec, pol_deg, pol_ang = self.model.interpolate(*params, integral_flux)
        energy = numpy.linspace(1., 12., 1000)
        for phase in numpy.linspace(0., 1., 100):
            min_pol_deg = pol_deg(energy, phase).min()
            msg = 'min pol deg = %.e at phase %.3f' % (min_pol_deg, phase)
            self.assertTrue(min_pol_deg >= 0., msg)

    def test_pol_ang_oscillations(self):
        """And, since the underlying models are noisy, we also have the
        issue of the angle swinging between -90 and 90---which is not only
        not exactly nice-looking, but also problematic when interpolating with
        a spline.
        """
        if xMagnetarTableModelT2020.auxfiles_missing():
            return
        params = (xMagnetarModelsT2020.SOLID_SURFACE_FIXED_IONS, 89.9, 60., 0.5, 0.34)
        integral_flux = 3.5e-11
        _, _, pol_ang_spline = self.model.interpolated_splines(*params, integral_flux)
        plt.figure('Position-angle oscillations')
        pol_ang_spline.plot()

    def _test_qedoff(self, file_name, *params):
        """Test the models with the QED effects disabled.
        """
        if xMagnetarTableModelT2020QedOff.auxfiles_missing():
            return
        table_I0, table_I, table_Q, table_U = self.model_qedoff._nearest_data(*params)
        # Load the reference data from the test file.
        file_path = os.path.join(IXPEOBSSIM_TEST_DATA, file_name)
        phase, energy_lo, energy_hi, I, Q, U = parse_data_file(file_path, *params[1:])
        phase /= (2. * numpy.pi)
        energy = 0.5 * (energy_lo + energy_hi)
        # Since in the model table we trim the first and last energy point,
        # we need to manually do the same here before we compare things.
        energy = energy[1:-1]
        I = I[:,1:-1]
        Q = Q[:,1:-1]
        U = U[:,1:-1]
        self.assertTrue(numpy.allclose(I, table_I0))
        self.assertTrue(numpy.allclose(Q, table_Q))
        self.assertTrue(numpy.allclose(U, table_U))
        spec_spline, pol_deg_spline, pol_ang_spline = self.model_qedoff.nearest_splines(*params, None)
        index = 0
        plt.figure('%s I' % file_name)
        plt.plot(energy, I[index], 'o')
        spec_spline.hslice(phase[index]).plot()
        setup_gca(logy=True)
        plt.figure('%s polarization degree' % file_name)
        pol_deg_spline.hslice(phase[index]).plot()
        plt.figure('%s polarization angle' % file_name)
        pol_ang_spline.hslice(phase[index]).plot()

    def test_qedoff_first(self):
        """
        """
        params = (xMagnetarModelsT2020.ATMOSPHERE, 0.1, 0.1, 0.3, 0.2)
        self._test_qedoff('atmosphere_0.10-0.10-03-02.dat', *params)

    def test_qedoff_random(self):
        """
        """
        params = (xMagnetarModelsT2020.SOLID_SURFACE_FREE_IONS, 179.9, 89.9, 1.4, 0.7)
        self._test_qedoff('solid_free_179.90-89.90-14-07.dat', *params)


    def test_print(self):
        """Basic printout test for low-level debugging.
        """
        if xMagnetarTableModelT2020.auxfiles_missing():
            return
        if xMagnetarTableModelT2020QedOff.auxfiles_missing():
            return
        file_path = os.path.join(IXPEOBSSIM_AUXFILES, xMagnetarTableModelT2020.I_FILE_NAME)
        hdu_list_qedon = fits.open(file_path)
        file_path = os.path.join(IXPEOBSSIM_AUXFILES, xMagnetarTableModelT2020QedOff.I_FILE_NAME)
        hdu_list_qedoff = fits.open(file_path)
        data_qedon = hdu_list_qedon['SPECTRA'].data
        data_qedoff = hdu_list_qedoff['SPECTRA'].data
        print('QED on parameter grid\n', data_qedon['PARAMVAL'])
        print('QED off parameter grid\n', data_qedoff['PARAMVAL'])
        self.assertTrue(numpy.allclose(data_qedon['PARAMVAL'], data_qedoff['PARAMVAL']))



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)

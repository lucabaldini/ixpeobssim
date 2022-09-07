#!/urs/bin/env python
#
# Copyright (C) 2020--2021, the ixpeobssim team.
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

"""Interface to the magnetar tabular models by Roberto Turolla and Roberto Taverna.
"""

from __future__ import print_function, division

import os
from enum import IntEnum

import numpy
from astropy.io import fits

from ixpeobssim import IXPEOBSSIM_AUXFILES, auxfiles_missing
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.core.fitsio import xPrimaryHDU
from ixpeobssim.core.fitting import power_law_analytical_fit, linear_analytical_fit
from ixpeobssim.core.spline import xInterpolatedBivariateSpline
from ixpeobssim.core.stokes import xModelStokesParameters
from ixpeobssim.utils.units_ import erg_to_keV, keV_to_erg
from ixpeobssim.utils.misc import pairwise

# pylint: disable=invalid-name, too-many-arguments, no-member, too-many-locals
# pylint: disable=too-many-instance-attributes, too-few-public-methods


def parse_legacy_data(file_path, pad_bounds=True, emin=1., emax=12.):
    """Parse a legacy data file, predating the glorious, full model table.

    .. warning::

        This is obsolete and, in principle, there should be no reason to use it
        except for testing.

    Here is an example of the underlying data format:

    .. code-block::

        89.9000      60.0000     0.500000     0.340000
        0   0.00100000     0.699132
        -0.300000    -0.273469      4456.25     0.776281      73.0709
        -0.273469    -0.246939      834.500     0.786651      73.3261

    The basic rules are:

    * the first line contains the model parameters: chi, xi, delta_phi, and beta;
    * then comes a variable number of blocks where the first (3-elements) row
      indicates the index and range of the bin phase and the following (5-element)
      rows encapsulate the spectral and polarimetric properties.

    More specifically, each 5-element row includes:

    * log10(energy_lo);
    * log10(energy_hi);
    * flux (arbitrary units---the number of counts in the underlying MC);
    * polarization degree;
    * polarization angle in decimal degrees.

    """
    phase, energy, flux, pol_deg, pol_ang = ([] for i in range(5))
    # Loop over the input file.
    with open(file_path, 'r') as input_file:
        input_file.readline()
        for line in input_file:
            values = [float(item) for item in line.split()]
            if len(values) == 3:
                _, lo, hi = values
                phase.append(0.5 * (lo + hi) / (2. * numpy.pi))
            if len(values) == 5:
                loglo, loghi, f, pd, pa = values
                e = 0.5 * (10.**loglo + 10.**loghi)
                if e < emin or e > emax:
                    continue
                energy.append(e)
                flux.append(f)
                pol_deg.append(pd)
                pol_ang.append(pa)
    # Convert the quantities to numpy arrays.
    num_phase_bins = len(phase)
    num_energy_bins = int(len(energy) / num_phase_bins)
    shape = (num_phase_bins, num_energy_bins)
    phase = numpy.array(phase)
    # Need only one copy of the array of energies...
    energy = numpy.array(energy[:num_energy_bins])
    # ...and reshape the others appropriately.
    flux = numpy.array(flux).reshape(shape)
    pol_deg = numpy.array(pol_deg).reshape(shape)
    pol_ang = numpy.array(pol_ang).reshape(shape)
    if pad_bounds:
        phase, flux, pol_deg, pol_ang = _pad_arrays(phase, flux, pol_deg, pol_ang)
    return phase, energy, flux, pol_deg, pol_ang



def _pad_arrays(phase, *arrays):
    """Pad the relevant arrays so that the resulting interpolator
    is well-behaved in 0 and 1.

    .. warning::

        This is only used for parsing legacy data and, in principle, there should
        be no reason to use it except for testing.
    """
    phase = numpy.hstack([phase[-1] - 1., phase, phase[0] + 1.])
    arrays = [numpy.vstack([array_[-1, :], array_, array_[0, :]]) for array_ in arrays]
    return [phase] + arrays



def parse_data_file(file_path, *target_params, num_phase_bins=9, num_energy_bins=49):
    """Parse the content of an input data file.

    The basic data structure is as follows:

    * a first line with the 4 parameter values;
    * a line with (phase_bin, phase_min, phase_max)
    * a series of lines with (logemin, logemax, I, Q/I, U/I)

    .. code-block::

             89.9000      89.9000      1.40000     0.200000
           1    0.0010000000      0.69913174
        -0.300000    -0.273469      852.250    0.0217493   -0.0900019
        -0.273469    -0.246939      725.750    0.0282958   -0.0814848
        ...
         0.946939     0.973469      78.0000    0.0389240   -0.0909807
         0.973469      1.00000      198.250    0.0437026   -0.0840274
           2      0.69913174       1.3972635
        ...
    """
    _parse_line = lambda line: tuple([float(item) for item in line.split()])
    logger.info('Parsing input file %s...', file_path)
    # Initialize the arrays for the output data.
    phase = numpy.zeros((num_phase_bins, ))
    energy_lo = numpy.zeros((num_energy_bins, ))
    energy_hi = numpy.zeros((num_energy_bins, ))
    I = numpy.zeros((num_phase_bins, num_energy_bins))
    Q = numpy.zeros((num_phase_bins, num_energy_bins))
    U = numpy.zeros((num_phase_bins, num_energy_bins))
    # Read the actual file.
    with open(file_path, 'r') as input_file:
        # Make sure the first lines agrees with the target parameters from the file name.
        line = input_file.readline()
        if len(target_params) > 0:
            assert _parse_line(line) == target_params
        phase_bin = 0
        energy_bin = num_energy_bins
        for line in input_file:
            # A line with 3 values, this signals the beginning of a new phase bin.
            values = _parse_line(line)
            if len(values) == 3:
                phase_index, phase_min, phase_max = values
                assert phase_index == phase_bin + 1
                assert energy_bin == num_energy_bins
                phase[phase_bin] = 0.5 * (phase_min + phase_max)
                phase_bin += 1
                energy_bin = 0
            # A line with 5 values contains the I, q and u values for a given energy bin.
            elif len(values) == 5:
                loge_min, loge_max, _I, _q, _u = values
                I[phase_bin - 1][energy_bin] = _I
                Q[phase_bin - 1][energy_bin] = _q * _I
                U[phase_bin - 1][energy_bin] = _u * _I
                if phase_bin == 1:
                    energy_lo[energy_bin] = 10.**loge_min
                    energy_hi[energy_bin] = 10.**loge_max
                else:
                    assert 10.**loge_min == energy_lo[energy_bin]
                    assert 10.**loge_max == energy_hi[energy_bin]
                energy_bin += 1
        assert phase_bin == num_phase_bins
    return phase, energy_lo, energy_hi, I, Q, U



def package_model_table(root_folder_path, model_name):
    """Build a model table package in OGIP-compliant FITS format given an archive
    of text files containing the Stokes spectra as a function of energy and
    phase.
    """
    logger.info('Packaging magnetar model table...')
    # Folder structure for the different physical models.
    folder_dict = {
        xMagnetarModelsT2020.ATMOSPHERE: 'atmosphere',
        xMagnetarModelsT2020.BLACKBODY: 'blackbody',
        xMagnetarModelsT2020.SOLID_SURFACE_FREE_IONS: 'solid(free)',
        xMagnetarModelsT2020.SOLID_SURFACE_FIXED_IONS: 'solid(fix)'
    }
    # Convenience function to construct the file name.
    def _file_name(chi, xi, delta_phi, beta):
        """Build the file name corresponding to a given set of parammeters.
        """
        delta_phi = '%02d' % (delta_phi * 10.)
        beta = ('%03d' % (beta * 100)).rstrip('0')
        return 'dataphaseen(%.2f-%.2f-%s-%s).dat' % (chi, xi, delta_phi, beta)

    # Loop over the input files and fill the spectra arrays.
    num_spectra = xParameterSpaceT2020.NUMBSPECTRA
    num_energy_bins = len(xParameterSpaceT2020.ENERGY_LO)
    num_phase_bins = len(xParameterSpaceT2020.PHASE)
    paramval = numpy.zeros((num_spectra, 6))
    intpspec_I = numpy.zeros((num_spectra, num_energy_bins))
    intpspec_Q = numpy.zeros((num_spectra, num_energy_bins))
    intpspec_U = numpy.zeros((num_spectra, num_energy_bins))
    logger.info('Looping over input files...')
    i = 0
    for model, folder_name in folder_dict.items():
        folder_path = os.path.join(root_folder_path, folder_name)
        for chi in xParameterSpaceT2020.CHI:
            for xi in xParameterSpaceT2020.XI:
                for delta_phi in xParameterSpaceT2020.DELTA_PHI:
                    for beta in xParameterSpaceT2020.BETA:
                        params = chi, xi, delta_phi, beta
                        file_name = _file_name(*params)
                        file_path = os.path.join(folder_path, file_name)
                        assert os.path.exists(file_path)
                        _phase, _elo, _ehi, I, Q, U = parse_data_file(file_path, *params)
                        assert numpy.allclose(_phase, xParameterSpaceT2020.PHASE)
                        assert numpy.allclose(_elo, xParameterSpaceT2020.ENERGY_LO)
                        assert numpy.allclose(_ehi, xParameterSpaceT2020.ENERGY_HI)
                        cols = (
                            numpy.full(num_phase_bins, model),
                            numpy.full(num_phase_bins, chi),
                            numpy.full(num_phase_bins, xi),
                            numpy.full(num_phase_bins, delta_phi),
                            numpy.full(num_phase_bins, beta),
                            _phase
                        )
                        _paramval = numpy.column_stack(cols)
                        idx = slice(i * num_phase_bins, (i + 1) * num_phase_bins)
                        paramval[idx,:] = _paramval
                        intpspec_I[idx,:] = I
                        intpspec_Q[idx,:] = Q
                        intpspec_U[idx,:] = U
                        i += 1
    logger.info('Done, %d file(s) parsed.', i)
    logger.info('Writing output files...')

    def _model_keywords(hduclas2=None, hduvers='1.1.0'):
        """Small nested function for the common keywords for the model headers.
        """
        keywords = [
            ('HDUCLASS', 'OGIP'),
            ('HDUCLAS1', 'XSPEC TABLE MODEL'),
            ('HDUVERS', hduvers),
        ]
        if hduclas2 is not None:
            keywords.insert(-2, ('HDUCLAS2', hduclas2))
        return keywords

    def _primary_hdu(model_name, model_type):
        """Small nested function for the primary header.

        See https://bitbucket.org/ixpesw/ixpeobssim/issues/438
        """
        header_keywords = _model_keywords() +\
            [
            ('MODLNAME', '%s_%s' % (model_name, model_type.lower()), 'Table model name'),
            ('MODLUNIT', '', 'Table model units'),
            ('REDSHIFT', False, 'Add redshift parameter?'),
            ('ADDMODEL', True, 'Is model additive?'),
            ('LOELIMIT', 0., 'Value of model below tabulated energies'),
            ('HIELIMIT', 0., 'Value of model above tabulated energies')
            ]
        return xPrimaryHDU(keywords=header_keywords)

    # PARAMETERS
    columns = [
        fits.Column(name='NAME', array=xParameterSpaceT2020.NAME, format='12A'),
        fits.Column(name='METHOD', array=numpy.zeros((6,), dtype=int), format='J'),
        fits.Column(name='INITIAL', array=xParameterSpaceT2020.INITIAL, format='E'),
        fits.Column(name='DELTA', array=xParameterSpaceT2020.DELTA, format='E'),
        fits.Column(name='MINIMUM', array=xParameterSpaceT2020.MINIMUM, format='E'),
        fits.Column(name='BOTTOM', array=xParameterSpaceT2020.MINIMUM, format='E'),
        fits.Column(name='TOP', array=xParameterSpaceT2020.MAXIMUM, format='E'),
        fits.Column(name='MAXIMUM', array=xParameterSpaceT2020.MAXIMUM, format='E'),
        fits.Column(name='NUMBVALS', array=xParameterSpaceT2020.NUMBVALS, format='J'),
        fits.Column(name='VALUE', array=xParameterSpaceT2020.VALUE, format='PE')
    ]
    parameters = fits.BinTableHDU.from_columns(columns, name='PARAMETERS')
    for key, value in _model_keywords('PARAMETERS', '1.0.0'):
        parameters.header.set(key, value)
    parameters.header.set('NINTPARM', 6)
    parameters.header.set('NADDPARM', 0)

    # ENERGIES
    columns = [
        fits.Column(name='ENERG_LO', array=xParameterSpaceT2020.ENERGY_LO, format='D'),
        fits.Column(name='ENERG_HI', array=xParameterSpaceT2020.ENERGY_HI, format='D')
    ]
    energies = fits.BinTableHDU.from_columns(columns, name='ENERGIES')
    for key, value in _model_keywords('ENERGIES', '1.0.0'):
        energies.header.set(key, value)

    # SPECTRA
    for spec, label in zip((intpspec_I, intpspec_Q, intpspec_U), ('I', 'Q', 'U')):
        columns = [
            fits.Column(name='PARAMVAL', array=paramval, format='6E'),
            fits.Column(name='INTPSPEC', array=spec, format='49E')
            ]
        spectra = fits.BinTableHDU.from_columns(columns, name='SPECTRA')
        for key, value in _model_keywords('SPECTRA', '1.0.0'):
            spectra.header.set(key, value)

        # Build the actual HDU list.
        args = model_name, label
        hdu_list = fits.HDUList([_primary_hdu(*args), parameters, energies, spectra])
        file_path = 'magnetar_rcs_model_table_2020MNRAS4925057T_qedoff_%s.fits' % label
        logger.info('Writing %s to %s...', label, file_path)
        hdu_list.writeto(file_path, overwrite=True)
    logger.info('Done.')



class xMagnetarModelsT2020(IntEnum):

    """Enum class encapsulating the physical models in Taverna et al. 2020
    https://ui.adsabs.harvard.edu/link_gateway/2020MNRAS.492.5057T/EPRINT_PDF

    * ``ATMOSPHERE = 1``
    * ``BLACKBODY = 2``
    * ``SOLID_SURFACE_FREE_IONS = 3``
    * ``SOLID_SURFACE_FIXED_IONS = 4``
    """

    ATMOSPHERE = 1
    BLACKBODY = 2
    SOLID_SURFACE_FREE_IONS = 3
    SOLID_SURFACE_FIXED_IONS = 4

    @classmethod
    def has_value(cls, value):
        """Convenience method that can be used downstrem to check the validity
        of a model passed as an integer.

        https://stackoverflow.com/questions/43634618
        """
        return value in cls._value2member_map_



class xParameterSpaceT2020:

    """Small namespace to encapsulate the standard grid of model parameters.
    """

    NAME = numpy.array(['Model', 'Chi (deg)', 'xi (deg)', 'delta_phi (r', 'beta', 'phase (rad)'])
    MODEL = [
        xMagnetarModelsT2020.ATMOSPHERE,
        xMagnetarModelsT2020.BLACKBODY,
        xMagnetarModelsT2020.SOLID_SURFACE_FREE_IONS,
        xMagnetarModelsT2020.SOLID_SURFACE_FIXED_IONS
    ]
    CHI = numpy.array(
        [0.1, 15, 30, 45, 60, 75, 89.9, 105, 120, 135, 150, 165, 179.9]
    )
    XI = numpy.array(
        [0.1, 15, 30, 45, 60, 75, 85, 89.9]
    )
    DELTA_PHI = numpy.array(
        [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
    )
    BETA = numpy.array(
        [0.2, 0.3, 0.34, 0.4, 0.5, 0.6, 0.7]
    )
    PHASE = numpy.array(
        [0.35006587, 1.04819762, 1.74632935, 2.4444611, 3.1425929, 3.84072455,
         4.53885635, 5.23698815, 5.9351197]
    )
    INITIAL = numpy.array([1., 60., 30., 0.5, 0.34, numpy.pi])
    DELTA = numpy.array([1., 5., 5., 0.1, 0.04, 0.1])
    MINIMUM = [min(item) for item in (MODEL, CHI, XI, DELTA_PHI, BETA, PHASE)]
    MAXIMUM = [max(item) for item in (MODEL, CHI, XI, DELTA_PHI, BETA, PHASE)]
    NUMBVALS = [len(item) for item in (MODEL, CHI, XI, DELTA_PHI, BETA, PHASE)]
    NUMBSPECTRA = numpy.prod(NUMBVALS)
    VALUE = [MODEL, CHI, XI, DELTA_PHI, BETA, PHASE]
    ENERGY_LO = numpy.array(
        [0.50118723, 0.53275925, 0.56631883, 0.60199377, 0.63991457, 0.68022564,
         0.72307609, 0.76862410, 0.81704297, 0.86851135, 0.92322190, 0.98137886,
         1.04319933, 1.10891409, 1.17876844, 1.25302345, 1.33195637, 1.41585898,
         1.50505025, 1.59986007, 1.70063847, 1.80776927, 1.92164429, 2.04269723,
         2.17137583, 2.30815516, 2.45355614, 2.60811057, 2.77240708, 2.94705336,
         3.13269418, 3.33003655, 3.53981039, 3.76279016, 3.99982509, 4.25178211,
         4.51962081, 4.80433187, 5.10696639, 5.42867697, 5.77064024, 6.13415859,
         6.52057658, 6.93132081, 7.36795560, 7.83207793, 8.32545543, 8.84991297,
         9.40738678]
    )
    ENERGY_HI = numpy.append(ENERGY_LO[1:], 10.)



class xMagnetarTableModelComponentT2020:

    """Interface to the magnetar tabular models.

    The models are coming as three separate FITS files for the I, Q, and U
    Stokes parameters, the basic structure of each one being, e.g.:

    .. code-block::

        No.    Name      Ver    Type      Cards   Dimensions   Format
          0  PRIMARY       1 PrimaryHDU      15   ()
          1  PARAMETERS    1 BinTableHDU     47   6R x 10C   [12A, J, E, E, E, E, E, E, J, PE(13)]
          2  ENERGIES      1 BinTableHDU     21   49R x 2C   [D, D]
          3  SPECTRA       1 BinTableHDU     21   314496R x 2C   [6E, 49E]

    (The actual dimensions of the cards might change, but the very fundamental
    structure of the files is fixed by the OGIP standard, I believe.)

    The ``PARAMETERS`` binary table encapsulate the following model parameters:

    * `Model`: the underlying Physical model, see the MagnetarModel enum
    * `chi`: angle between the line of sight and the rotation angle of the star.
    * `xi`: angle between the magnetic axis ans the rotation angle of the star.
    * `delta_phi`: twist angle.
    * `beta`: electron velocity in units of c.
    * `phase`: the phase bins.

    For illustration purposes, in the initial implementation we have
    4 physical models x 13 chi x 8 xi x 12 delta_phi x 7 beta x 9 phase bins,
    that is, 314496 combinations. Spelling things out explicitely:

    .. code-block::

        Model     -> [1, 2, 3, 4]
        chi       -> [0.1, 15, 30, 45, 60, 75, 89.9, 105, 120, 135, 150, 165, 179.9]
        xi        -> [0.1, 15, 30, 45, 60, 75, 85, 89.9]
        delta_phi -> [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
        beta      -> [0.2, 0.3, 0.34, 0.4, 0.5, 0.6, 0.7]
        phase     -> [0.350066, 1.0482, 1.74633, 2.44446, 3.14259, 3.84072,
                      4.53886, 5.23699, 5.93512]

    The ``ENERGY`` binary extensions encapsulates the energy binning, and in the
    initial implementations has 49 bins between ~0.5 and 10 keV.

    Finally, the glorious ``SPECTRA`` extensions contains the spectra for all the
    possible configurations. It comes with two columns, the first of which
    contains the 6 model parameters and the second the 49 spectral values.

    Note the ``SPECTRA`` extension is filled by looping over the model parameters in
    reverse order, i.e., from `phase` to `Model`.
    """

    def __init__(self, file_path, normalize=False):
        """Open the tabular model FITS files and cache the necessary data for use.
        """
        logger.info('Loading tabular models from %s...', file_path)
        with fits.open(file_path) as hdu_list:
            hdu_list.info()
            # Cache the necessary stuff so that we can close the underlying
            # FITS file and move on.
            self.par_names = tuple(hdu_list['PARAMETERS'].data['NAME'])
            self.par_values = tuple(hdu_list['PARAMETERS'].data['VALUE'])
            self.energ_lo = hdu_list['ENERGIES'].data['ENERG_LO']
            self.energ_hi = hdu_list['ENERGIES'].data['ENERG_HI']
            self.spectral_params = hdu_list['SPECTRA'].data['PARAMVAL']
            self.spectral_data = hdu_list['SPECTRA'].data['INTPSPEC']
        # Normalize by the bin width in energy---this is needed for the I
        # component.
        if normalize:
            self.spectral_data /= (self.energ_hi - self.energ_lo)
        # Calculate the shape of the table---this is essentially the number
        # of values each of the parameters can take, e.g., (4, 13, 8, 12, 7, 9).
        self.shape = tuple([len(_value) for _value in self.par_values])
        # Auxiliary variable to retrieve the spectra by parameter index.
        # This is essentially the cumulative product of the last part of the
        # shape in reverse order, e.g., (78624, 6048, 756, 63, 9).
        self._cumshape = numpy.flip(numpy.cumprod(self.shape[::-1]))[1:]
        # Is this the right thing? Need to talk with Roberto^2 and see
        # what they are actually doing when they generate the table.
        self.energy = 0.5 * (self.energ_lo + self.energ_hi)
        # And, sice we're at it, cache the phase grid for later use. Note that
        # the phase in the model runs from 0 to 2pi, while for the ixpeobsim
        # purpopses we do need a number between 0 and 1.
        self.phase = self.par_values[-1] / (2. * numpy.pi)
        # Trim the first and the last energy points in the energy grid and all
        # the related spectra, as in the underlying Monte Carlo simulation
        # the outesrmost bins act as underflow/overflow accumulators, and
        # the corresponding values are therefore unphysical.
        self.energy = self.energy[1:][:-1]
        self.spectral_data = self.spectral_data[:, 1:][:, :-1]
        logger.info('Done.\n%s', self)

    def __str__(self):
        """String formatting.
        """
        text = 'Table shape %s\n' % str(self.shape)
        for name, values in zip(self.par_names, self.par_values):
            text += '%s: %s\n' % (name, values)
        return text.strip('\n')

    def num_phase_bins(self):
        """Return the number of phase bins.
        """
        return self.shape[-1]

    def parameters(self, *indices):
        """Retrieve the parameter values for a given set of indices.
        """
        return tuple(self.par_values[i][j] for i, j in enumerate(indices))

    def spectrum_data(self, *indices):
        """Retrieve the spectrum data for a given series of indices.

        Note that the function requires to pass *all* the indices, *except* for
        the last one, i.e., the one referring to the bin phase, and the
        return value contains the data for all the bin phases.
        """
        assert len(indices) == len(self.shape) - 1
        index = (numpy.array(indices) * self._cumshape).sum()
        slice_ = slice(index, index + self.num_phase_bins())
        return self.spectral_data[slice_]

    @staticmethod
    def _nearest_index(array_, value):
        """Find the index of the nearest array element to the target value.

        Note that, since we are only using this to bisect the arrays of the
        parameter values, we're not doing a binary search. See
        https://stackoverflow.com/questions/2566412
        for a more extensive discussion.
        """
        return numpy.abs(array_ - value).argmin()

    def nearest_indices(self, *params):
        """Return the indices of the slice of the underlying multi-dimensional
        array corresponding to the model that comes closest to the input parameters.
        """
        par_values = self.par_values[:-1]
        return tuple([self._nearest_index(a, v) for a, v in zip(par_values, params)])

    def nearest(self, *params):
        """Return the slice of the underlying multi-dimensional array
        corresponding to the model that come closest to the input parameters.
        """
        indices = self.nearest_indices(*params)
        return self.spectrum_data(*indices)

    @staticmethod
    def _combinations(values):
        """Build an array with all the combinations of a given set of iterables.

        This is used when interpolating the model table, where we bisect the
        all the parameter grids and we have to build a weighted average of all
        possible combinations.

        See https://stackoverflow.com/questions/1208118 for more details.
        """
        return numpy.array(numpy.meshgrid(*values)).T.reshape(-1, len(values))

    def interpolate_indices(self, *params, threshold=1.e-6):
        """Interpolate the underlying SPECTRA extension at a given set of input
        parameters.

        This performs a simple linear interpolation, where the first parameter
        (the actual model) is treated in a special fashion in that is not
        interpolated at all---in fact the function requires that the value
        passed as an argument to the function is one of the four valis values.

        For all the other parameters, the corresponding grids of tabulated values
        are bisected at the target value and the two nearest elements are
        weighted proportinally to the relative distance of the target.
        """
        model, *params = params
        model = int(model)
        assert xMagnetarModelsT2020.has_value(model)
        par_values = self.par_values[1:-1]
        indices = []
        weights = []
        # Bisect the parameter grids.
        for value, grid in zip(params, par_values):
            i = numpy.searchsorted(grid, value)
            indices.append((i - 1, i))
            lo, hi = grid[i - 1], grid[i]
            w = (hi - value) / (hi - lo)
            # Small hack to avoid numerical artifacts due to the fact that the
            # parameter grids are saved in single precision, while everywhere
            # else we are operating in double precision. Essentially if the
            # weights are too close to either 0 or 1, we attribute that to simple
            # numerical rounding.
            if w < threshold:
                w = 0.
            elif 1. - w < threshold:
                w = 1.
            weights.append((w, 1. - w))
        # Create all the combination of the relevant parameters and weights.
        indices = self._combinations(indices)
        weights = self._combinations(weights)
        # Multiply the weights along the proper axis to come up with one value
        # for each parameter set.
        weights = numpy.prod(weights, 1)
        # Get rid of the combinations that will not contribute to the weighted
        # average.
        mask = weights > 0.
        indices = indices[mask]
        weights = weights[mask]
        # Prepend the model to the indices.
        model_indices = numpy.full((len(indices), 1), model - 1)
        indices = numpy.hstack((model_indices, indices))
        return indices, weights

    def weighted_average(self, indices, weights):
        """Calculate the weighted average for a given set of entries in the
        underlying SPECTRA table.
        """
        data = numpy.full((self.shape[-1], len(self.energy)), 0.)
        for i, w in zip(indices, weights):
            data += self.spectrum_data(*i) * w
        return data

    def interpolate(self, *params):
        """Interpolate the underlying SPECTRA table to a given set of parameters.
        """
        indices, weights = self.interpolate_indices(*params)
        return self.weighted_average(indices, weights)



class xMagnetarTableModelT2020:

    """Full description of a magnetar spectro-polarimetric model.
    """

    I_FILE_NAME = 'magnetar_rcs_model_table_2020MNRAS4925057T_I.fits'
    Q_FILE_NAME = 'magnetar_rcs_model_table_2020MNRAS4925057T_Q.fits'
    U_FILE_NAME = 'magnetar_rcs_model_table_2020MNRAS4925057T_U.fits'

    def __init__(self):
        """Constructor.
        """
        if self.auxfiles_missing():
            abort()
        # Note we need to keep track of the I spectrum in two flavors:
        # * I0: the un-normalized spectrum, straight from the input tables,
        #   that is needed to normalize Q and U downstream;
        # * I: the normalized spectrum, where the bin contents are divided by
        #   the bin width to get the right spectrum in cm^{-2} s^{-1} keV^{-1}
        #
        # See  https://bitbucket.org/ixpesw/ixpeobssim/issues/453
        # for more details about the matter.
        self._I0 = self._load_component(self.I_FILE_NAME)
        self._I = self._load_component(self.I_FILE_NAME, normalize=True)
        self._Q = self._load_component(self.Q_FILE_NAME)
        self._U = self._load_component(self.U_FILE_NAME)
        # Minimal check on the compatibility of the components.
        assert numpy.allclose(self._I.energy, self._Q.energy)
        assert numpy.allclose(self._I.energy, self._U.energy)
        assert numpy.allclose(self._I.phase, self._Q.phase)
        assert numpy.allclose(self._I.phase, self._U.phase)
        # Cache the params for the energy and phase binning.
        self.energy = self._I.energy
        self.phase = self._I.phase
        self.par_names = self._I.par_names
        self.par_values = self._I.par_values

    @staticmethod
    def _load_component(file_name, normalize=False):
        """Load a model component.
        """
        file_path = os.path.join(IXPEOBSSIM_AUXFILES, file_name)
        return xMagnetarTableModelComponentT2020(file_path, normalize)

    @classmethod
    def auxfiles_missing(cls):
        """Convenience function to check for missing auxiliary files.

        Note this is a classmethod, so it can be invoked without an instance
        of the class, which can be handy, e.g., for skipping unit tests.
        """
        return auxfiles_missing(cls.I_FILE_NAME, cls.Q_FILE_NAME, cls.U_FILE_NAME)

    def _nearest_data(self, model, chi, xi, delta_phi, beta):
        """Finde the underlying spectrum that comes closest to the input parameters.
        """
        assert xMagnetarModelsT2020.has_value(model)
        params = model, chi, xi, delta_phi, beta
        indices = self._I.nearest_indices(*params)
        return self._I0.spectrum_data(*indices),\
            self._I.spectrum_data(*indices),\
            self._Q.spectrum_data(*indices),\
            self._U.spectrum_data(*indices)

    def _interpolate_data(self, model, chi, xi, delta_phi, beta):
        """Interpolate the underlying Stokes spectra.
        """
        assert xMagnetarModelsT2020.has_value(model)
        params = model, chi, xi, delta_phi, beta
        indices, weights = self._I.interpolate_indices(*params)
        return self._I0.weighted_average(indices, weights),\
            self._I.weighted_average(indices, weights),\
            self._Q.weighted_average(indices, weights),\
            self._U.weighted_average(indices, weights)

    @staticmethod
    def _pad_phase_data(data):
        """Small conveinence function to pad spectral data in phase, see the
        documentation of xMagnetarTableModelT2020._pad_phase() for details.
        """
        data0 = 0.5 * (data[0] + data[-1])
        return numpy.concatenate(([data[-1], data0], data, [data0, data[0]]))

    @staticmethod
    def _pad_phase(phase, I0, I, Q, U):
        """Pad the phase grid and the Stokes data to guarantee a proper
        continuity condition at the bounduaries.

        Mond the underlying phase grid runs from the the center of the first
        phase bin (> 0.) to the center of the last phase bin (< 1.), so we have to
        do something to avoid bogus extrapolations near 0 and 1.
        The basic paddind strategy is to add two points on the left, at
        max(phase) - 1. and 0., and two points to the right, at 1. and min(phase) + 1.,
        so that the new grids covers the range max(phase) - 1. (< 0) to
        min(phase) + 1. (> 1.). The spectral data are padded accordingly, averaging
        the first and last phase points to set the value at 0. and 1.

        Note that, with this prescription, using a spline of order k = 2 for the
        phase axis when building the relevant 2-dimensional spline guarantees that
        the values at 0. and 1. are identical---i.e., you have no jumps when
        the phase rolls over.
        """
        phase = numpy.concatenate(([phase.max() - 1., 0.], phase, [1., phase.min() + 1.]))
        I0 = xMagnetarTableModelT2020._pad_phase_data(I0)
        I = xMagnetarTableModelT2020._pad_phase_data(I)
        Q = xMagnetarTableModelT2020._pad_phase_data(Q)
        U = xMagnetarTableModelT2020._pad_phase_data(U)
        return phase, I0, I, Q, U

    @staticmethod
    def _pad_energy(energy, I, q, u, emax=15., num_points=10, fit_emin=6., fit_emax=10.):
        """Pad the energy grid and Stokes data above the maximum energy (~9 keV)
        so we do not rely on the naive (constrant) 2-d spline extrapolation.

        This turned out to be more complicated than originally anticipates as,
        while it's reasonably easy to fit the high-energy part of I with a power law
        and use that to extrapolate, this doesn't really work for Q and U, that
        can go negative.

        The strategy we adopt is a two-step one, namely:

        * we fit I in each phase bin with a power law between fit_emin and
          fit_emax, and we use the fit to extrapolate above fit_emax;
        * we fit q = Q/I and u = U/I in each phase bin with a straight line
          between fit_emin and fit_emax, and we use the fit to extrapolate above
          fit_emax.
        """
        num_phase_bins, _ = I.shape
        mask = numpy.logical_and(energy >= fit_emin, energy <= fit_emax)
        pad_energy = numpy.linspace(fit_emax, emax, num_points)
        # Small anonymous function to select the part of the I, Q or U
        # two-dimensional arrays to be used for the fitting, which requires some
        # numpy magic.
        mask2d = numpy.tile(mask, (num_phase_bins, 1))
        _subselect = lambda data: data[mask2d].reshape((num_phase_bins, len(x)))

        # First step: fit the high-energy part of the I data with a power law
        # and pad them with the fit extrapolation.
        x = energy[mask]
        y = _subselect(I).T
        norm, index = power_law_analytical_fit(x, y)
        norm = numpy.tile(norm, (num_points, 1)).T
        index = numpy.tile(index, (num_points, 1)).T
        I_padding = norm * pad_energy**index

        # Second step: extrapolate Q/I and U/I.
        def _qupadding(data):
            """Small convenience function to fit Q and U and extrapolate to
            the lower-end of the energy padding.
            """
            y = _subselect(data).T
            m, q = linear_analytical_fit(x, y)
            m = numpy.tile(m, (num_points, 1)).T
            q = numpy.tile(q, (num_points, 1)).T
            _pad = m * pad_energy + q
            return _pad

        q = numpy.hstack((q, _qupadding(q)))
        u = numpy.hstack((u, _qupadding(u)))
        # Note that I goes last so that we have access to the original array
        # up to the last minute, and we can use it, e.g., to estimate a proxy
        # of the errors in _qupadding.
        I = numpy.hstack((I, I_padding))
        energy = numpy.hstack((energy, numpy.linspace(fit_emax, emax, num_points)))
        return energy, I, q, u

    @staticmethod
    def _build_splines(energy, phase, I0, I, Q, U, integral_flux, emin, emax,
                       pad_phase=True, pad_energy=True):
        """Build a set of three splines that can be used by ixpeobssim to run a
        simulation.

        Note that, With our prescrition on the padding, using a spline of order
        2 for the phase axis guarantees continuity when the phase rolls over,
        so we don't expose ky in the public API.

        Likewise, the spline order over the energy axis is hard-coded.

        The flags for disabling the padding in phase and energy are here mainly
        for testing purposes---so that we can disengage the padding and compare
        the spline interpolation with the underlying data directly.
        """
        kx = 2
        ky = 2
        if pad_phase:
            phase, I0, I, Q, U = xMagnetarTableModelT2020._pad_phase(phase, I0, I, Q, U)
        q = Q / I0
        u = U / I0
        if pad_energy:
            energy, I, q, u = xMagnetarTableModelT2020._pad_energy(energy, I, q, u)
        pd = xModelStokesParameters.polarization_degree(q, u)
        pa = xModelStokesParameters.polarization_angle(q, u)
        fmt = dict(kx=kx, ky=ky, xlabel='Energy [keV]', ylabel='Phase',
                   zlabel='Flux [cm$^{-2}$ s$^{-1}$ keV$^{-1}$]')
        spec_spline = xInterpolatedBivariateSpline(energy, phase, I.T, **fmt)
        if integral_flux is not None:
            # Normalize the spectral spline to the correct integral flux.
            # Note that we need a custom spline for the energy flux, here---see
            # https://bitbucket.org/ixpesw/ixpeobssim/issues/471
            espec_spline = xInterpolatedBivariateSpline(energy, phase, (energy * I).T)
            norm = erg_to_keV(integral_flux / espec_spline.integral(emin, emax, 0., 1.))
            spec_spline = spec_spline.scale(norm)
        fmt.update({'zlabel': 'Polarization degree'})
        pol_deg_spline = xInterpolatedBivariateSpline(energy, phase, pd.T, **fmt)
        fmt.update({'zlabel': 'Polarization angle [rad]'})
        pol_ang_spline = xInterpolatedBivariateSpline(energy, phase, pa.T, **fmt)
        return spec_spline, pol_deg_spline, pol_ang_spline

    @staticmethod
    def _wrap_splines(spec_spline, pol_deg_spline, pol_ang_spline):
        """Wrap a set of splines in such a way they have the right signature
        for ixpeobssim to use them.

        Note that we do slightly more than just wrapping the splines, as,
        for the case of the polarization degree, we also make sure that
        the values returned lie within the physical range [0., 1.]. This is
        necessary because the underlying models are noisy, and when
        interpolating with a cubic spline it is possible to undershoot or
        overshoot the physical bounds.
        """
        def spec(E, t):
            """Nested spectrum definition.
            """
            return spec_spline(E, t)

        def pol_deg(E, t, ra=None, dec=None):
            """Nested polarization degree definition.

            Note we are clipping the values to the [0., 1.] interval.
            """
            # pylint: disable=unused-argument
            return numpy.clip(pol_deg_spline(E, t), 0., 1.)

        def pol_ang(E, t, ra=None, dec=None):
            """Nested polarization angle definition.
            """
            # pylint: disable=unused-argument
            return pol_ang_spline(E, t)

        return spec, pol_deg, pol_ang

    def nearest_splines(self, model, chi, xi, delta_phi, beta, integral_flux,
        emin=2., emax=10., pad_phase=True, pad_energy=True):
        """Return a set of model splines for the nearest set of parameters in the
        underlying model table.
        """
        data = self._nearest_data(model, chi, xi, delta_phi, beta)
        args = self.energy, self.phase, *data, integral_flux, emin, emax, pad_phase, pad_energy
        return self._build_splines(*args)

    def interpolated_splines(self, model, chi, xi, delta_phi, beta, integral_flux,
        emin=2., emax=10., pad_phase=True, pad_energy=True):
        """Return a set of model splines interpolating the underlying model table
        to the target parammeters.
        """
        data = self._interpolate_data(model, chi, xi, delta_phi, beta)
        args = self.energy, self.phase, *data, integral_flux, emin, emax, pad_phase, pad_energy
        return self._build_splines(*args)

    def nearest(self, model, chi, xi, delta_phi, beta, integral_flux, emin=2., emax=10.):
        """Return a set of wraps of the nearest splines ready to be used in ixpeobssim.
        """
        args = model, chi, xi, delta_phi, beta, integral_flux, emin, emax
        splines = self.nearest_splines(*args)
        return self._wrap_splines(*splines)

    def interpolate(self, model, chi, xi, delta_phi, beta, integral_flux, emin=2., emax=10.):
        """Return a set of wraps of the interpolated splines ready to be used in ixpeobssim.
        """
        args = model, chi, xi, delta_phi, beta, integral_flux, emin, emax
        splines = self.interpolated_splines(*args)
        return self._wrap_splines(*splines)

    def energy_spectrum(self, model, chi, xi, delta_phi, beta, integral_flux, emin=2., emax=10.):
        """Return a bivariate spline representing the energy spectrum.

        This is essentially retrieving the photon spectrum (as a bivariate spline)
        and multiplying the z array with the proper tiling of the x array
        (i.e., the energy).
        """
        args = model, chi, xi, delta_phi, beta, integral_flux, emin, emax
        spec_spline, _, _ = self.interpolated_splines(*args)
        z = spec_spline.z * numpy.tile(spec_spline.x, (len(spec_spline.y), 1)).T
        return xInterpolatedBivariateSpline(spec_spline.x, spec_spline.y, z)

    def phase_resolved_integral_flux(self, phase_binning, model, chi, xi, delta_phi,
        beta, integral_flux, emin=2., emax=10.):
        """Return the integral flux of the model (in erg/cm^2) in arbitrary phase bins.
        """
        espec = self.energy_spectrum(model, chi, xi, delta_phi, beta, integral_flux, emin, emax)
        flux = numpy.array([espec.integral(emin, emax, min_, max_) / (max_ - min_) \
            for (min_, max_) in pairwise(phase_binning)])
        return keV_to_erg(flux)



class xMagnetarTableModelT2020QedOff(xMagnetarTableModelT2020):

    """Full description of a magnetar spectro-polarimetric model no QED effects.
    """

    I_FILE_NAME = 'magnetar_rcs_model_table_2020MNRAS4925057T_qedoff_I.fits'
    Q_FILE_NAME = 'magnetar_rcs_model_table_2020MNRAS4925057T_qedoff_Q.fits'
    U_FILE_NAME = 'magnetar_rcs_model_table_2020MNRAS4925057T_qedoff_U.fits'




if __name__ == '__main__':
    package_model_table('/data/temp/qedoff', 'Magnetar_1708_qedoff')

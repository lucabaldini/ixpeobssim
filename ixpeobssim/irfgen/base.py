#!/usr/bin/env python
#
# Copyright (C) 2021--2022, the ixpeobssim team.
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

"""Base module for all the facilities related to the creation of response functions.
"""


from __future__ import print_function, division

from astropy.io import fits
import numpy

from ixpeobssim.core.fitsio import xPrimaryHDU
from ixpeobssim.instrument.du import du_physical_name, du_logical_name
from ixpeobssim.irf.caldb import irf_file_path, irf_file_name, parse_irf_name
from ixpeobssim.irf.ebounds import ENERGY_GRID, ENERGY_LO, ENERGY_HI
from ixpeobssim.irf.ebounds import NUM_CHANNELS, TLMIN, TLMAX, channel_to_energy
from ixpeobssim.irfgen.fmt import xBinTableHDUSPECRESPARF, xBinTableHDUSPECRESPMRF,\
    xBinTableHDUSPECRESPMODF, xBinTableHDUVIGNETTING, xBinTableHDUPSF, \
    xBinTableHDUMATRIX, xBinTableHDUEBOUNDS
from ixpeobssim.irfgen.gpd import quantum_efficiency, modulation_factor, xEdispDataInterface
from ixpeobssim.irfgen.gpd import GPD_FILL_TEMPERATURE
from ixpeobssim.irfgen.mma import effective_area, vignetting_spline
from ixpeobssim.irfgen.du import uv_filter_transparency
from ixpeobssim.utils.logging_ import logger


# pylint: disable=invalid-name, too-many-arguments, too-many-locals

# Version of the auxiliary data files to be used for the IRF creation.
AUX_VERSION = 3


def _primary_hdu(irf_type, du_id, creator, comments):
    """Small wrapper function for creating the primary HDU.
    """
    keywords = [
        ('TELESCOP', 'IXPE', 'Telescope name'),
        ('INSTRUME', 'GPD', 'Instrument name'),
        ('DETNAM'  , du_logical_name(du_id), 'Name of the logical detector unit'),
        ('DET_ID'  , du_physical_name(du_id), 'Name of the physical detector unit'),
        ('IRFTYPE', irf_type, 'type of the file')
    ]
    hdu = xPrimaryHDU(creator=creator, keywords=keywords, comments=comments)
    logger.info('Creating PRIMARY HDU...')
    return hdu, keywords


def _update_bin_table_header(bin_table, irf_name, irf_type, du_id, weight_name, ccnm0001, cdes0001):
    """Update the binary table extension header with all the keyword
    values that cannot be inferred at creation time.

    This was added to align the file format to the IXPE CALDB, in response to
    https://bitbucket.org/ixpesw/ixpeobssim/issues/468
    """
    base, intent, version = parse_irf_name(irf_name)
    file_name = irf_file_name(base, du_id, irf_type, intent, version)
    weight_name = str(weight_name).upper()
    bin_table.header.set('VERSION', version)
    bin_table.header.set('FILENAME', file_name)
    bin_table.header.set('CCNM0001', ccnm0001)
    bin_table.header.set('CDES0001', cdes0001)
    bin_table.header.set('CBD10001', 'WEIGHT(%s)' % weight_name)


def _write_irf(hdu_list, irf_name, du_id, irf_type):
    """Small wrapper function to write an HDU list to file.
    """
    file_path = irf_file_path(irf_name, du_id, irf_type, check_file=False)
    logger.info('Writing %s data to %s...', irf_type, file_path)
    hdu_list.info()
    hdu_list.writeto(file_path, overwrite=True, checksum=True)
    hdu_list.close()
    logger.info('Done.')


def _aeff_data(du_id, pressure, weight_name=None, aux_version=AUX_VERSION):
    """Retrieve the effective area data.

    This is factored out, as we need the information both for the arf and for the
    mrf files.
    """
    aeff = effective_area(ENERGY_GRID, du_id)
    aeff *= uv_filter_transparency(ENERGY_GRID)
    args = pressure, weight_name, aux_version
    aeff *= quantum_efficiency(ENERGY_GRID, GPD_FILL_TEMPERATURE, *args)
    return aeff


def make_arf(irf_name, du_id, pressure, creator, comments, weight_name=None,
    aux_version=AUX_VERSION):
    """Write the IXPE effective area response function.
    """
    irf_type = 'arf'
    logger.info('Creating %s FITS file for DU %s...', irf_type, du_id)
    aeff = _aeff_data(du_id, pressure, weight_name, aux_version)
    primary_hdu, header_keywords = _primary_hdu(irf_type, du_id, creator, comments)
    data = [ENERGY_LO, ENERGY_HI, aeff]
    specresp_hdu = xBinTableHDUSPECRESPARF(data, header_keywords, comments)
    ccnm0001 = 'SPECRESP'
    cdes0001 = 'IXPE DU%d Ancillary Response File' % du_id
    _update_bin_table_header(specresp_hdu, irf_name, irf_type, du_id, weight_name, ccnm0001, cdes0001)
    hdu_list = fits.HDUList([primary_hdu, specresp_hdu])
    _write_irf(hdu_list, irf_name, du_id, irf_type)


def make_mrf(irf_name, du_id, pressure, creator, comments, weight_name=None,
    aux_version=AUX_VERSION):
    """Write the IXPE effective area response function.
    """
    irf_type = 'mrf'
    logger.info('Creating %s FITS file for DU %s...', irf_type, du_id)
    mrf = _aeff_data(du_id, pressure, weight_name, aux_version)
    mrf *= modulation_factor(ENERGY_GRID, pressure, weight_name, aux_version)
    primary_hdu, header_keywords = _primary_hdu(irf_type, du_id, creator, comments)
    data = [ENERGY_LO, ENERGY_HI, mrf]
    specresp_hdu = xBinTableHDUSPECRESPMRF(data, header_keywords, comments)
    ccnm0001 = 'MODSPECRESP'
    cdes0001 = 'IXPE DU%d Modulation Response Function' % du_id
    _update_bin_table_header(specresp_hdu, irf_name, irf_type, du_id, weight_name, ccnm0001, cdes0001)
    hdu_list = fits.HDUList([primary_hdu, specresp_hdu])
    _write_irf(hdu_list, irf_name, du_id, irf_type)


def make_vign(irf_name, du_id, pressure, creator, comments, weight_name=None,
    aux_version=AUX_VERSION):
    """Write the IXPE vignetting response function.
    """
    #pylint: disable=unused-argument
    irf_type = 'vign'
    logger.info('Creating %s FITS file for DU %s...', irf_type, du_id)
    vign_spline = vignetting_spline()
    theta = vign_spline.y
    vignetting = vign_spline(*numpy.meshgrid(ENERGY_GRID, theta))
    primary_hdu, header_keywords = _primary_hdu(irf_type, du_id, creator, comments)
    data = [ENERGY_LO, ENERGY_HI, theta, vignetting]
    vign_hdu = xBinTableHDUVIGNETTING(data, header_keywords, comments)
    hdu_list = fits.HDUList([primary_hdu, vign_hdu])
    _write_irf(hdu_list, irf_name, du_id, irf_type)


def make_rmf(irf_name, du_id, pressure, creator, comments, weight_name=None,
    aux_version=AUX_VERSION):
    """Write the IXPE edisp response function.

    The specifications are described at page ~15 of the following document:
    ftp://legacy.gsfc.nasa.gov/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.ps
    """
    irf_type = 'rmf'
    logger.info('Creating %s FITS file for DU %s...', irf_type, du_id)
    primary_hdu, header_keywords = _primary_hdu(irf_type, du_id, creator, comments)
    # Update the header keywords with the bits that are necessary for the rmf.
    header_keywords += \
        [('DETCHANS', NUM_CHANNELS, 'Total number of detector channels'),
         ('TLMIN4', TLMIN, 'First channel number'),
         ('TLMAX4', TLMAX, 'Last channel number')]
    num_rows = len(ENERGY_LO)
    ngrp = numpy.ones(num_rows)
    fchan = numpy.zeros(num_rows)
    nchan = numpy.array([NUM_CHANNELS] * num_rows, 'i')
    edisp = xEdispDataInterface(weight_name, aux_version)
    hist = edisp.combine(GPD_FILL_TEMPERATURE, pressure)
    matrix = hist.content.T
    assert (matrix >= 0.).all()
    ch = numpy.arange(NUM_CHANNELS)
    data = [ENERGY_LO, ENERGY_HI, ngrp, fchan, nchan, matrix]
    matrix_hdu = xBinTableHDUMATRIX(NUM_CHANNELS, data, header_keywords, comments)
    ccnm0001 = 'MATRIX'
    cdes0001 = 'IXPE DU%d Response Matrix' % du_id
    _update_bin_table_header(matrix_hdu, irf_name, irf_type, du_id, weight_name, ccnm0001, cdes0001)
    ch = numpy.arange(NUM_CHANNELS)
    emin = channel_to_energy(ch - 0.5)
    emax = channel_to_energy(ch + 0.5)
    data = [ch, emin, emax]
    ebounds_hdu = xBinTableHDUEBOUNDS(data, header_keywords, comments)
    ccnm0001 = 'EBOUNDS'
    cdes0001 = 'IXPE DU%d Response Matrix' % du_id
    _update_bin_table_header(ebounds_hdu, irf_name, irf_type, du_id, weight_name, ccnm0001, cdes0001)
    hdu_list = fits.HDUList([primary_hdu, matrix_hdu, ebounds_hdu])
    _write_irf(hdu_list, irf_name, du_id, irf_type)


def make_modf(irf_name, du_id, pressure, creator, comments, weight_name=None,
    aux_version=AUX_VERSION):
    """Write the IXPE modulation factor response function.
    """
    irf_type = 'modf'
    logger.info('Creating %s FITS file for DU %s...', irf_type, du_id)
    modf = modulation_factor(ENERGY_GRID, pressure, weight_name, aux_version)
    primary_hdu, header_keywords = _primary_hdu(irf_type, du_id, creator, comments)
    data = [ENERGY_LO, ENERGY_HI, modf]
    modf_hdu = xBinTableHDUSPECRESPMODF(data, header_keywords, comments)
    ccnm0001 = 'MODFACT'
    cdes0001 = 'IXPE DU%d Modulation Factor' % du_id
    _update_bin_table_header(modf_hdu, irf_name, irf_type, du_id, weight_name, ccnm0001, cdes0001)
    hdu_list = fits.HDUList([primary_hdu, modf_hdu])
    _write_irf(hdu_list, irf_name, du_id, irf_type)


def make_psf(irf_name, du_id, params, creator, comments, weight_name=None,
    aux_version=AUX_VERSION):
    """Write the IXPE PSF response function.
    """
    #pylint: disable=unused-argument
    irf_type = 'psf'
    logger.info('Creating %s FITS file for DU %s...', irf_type, du_id)
    params = [numpy.array([par]) for par in params]
    primary_hdu, header_keywords = _primary_hdu(irf_type, du_id, creator, comments)
    psf_hdu = xBinTableHDUPSF(params, header_keywords, comments)
    hdu_list = fits.HDUList([primary_hdu, psf_hdu])
    _write_irf(hdu_list, irf_name, du_id, irf_type)

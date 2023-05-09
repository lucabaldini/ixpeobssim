#!/usr/bin/env python
#
# Copyright (C) 2021--2023, the ixpeobssim team.
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

__description__ = \
"""Format reconstructed event lists generated with ixpesim to make them
interoperable with ixpeobssim (e.g., for spectro-polarimetric fits in XSPEC).
Phrased in a slightly different way, this is tranforming the rough equivalent
of Level-1 files produced by ixpesim into the rough equivalent of Level-2 files.

This is adding all the necessary keywords to the relevant headers, as well as
a few columns in the EVENTS extensions that are necessary for the
spectro-polarimetric analysis in XSPEC, most notably:
- PI
- Q
- U
- RA, DEC
- X, Y
- W_MOM

Additionally, the script allows to strip off useless columns, producing smaller
output files that easier to exchange in several different context. This
functionality is controlled by the --stripmode command-line switch and, more
specifically:

* --stripmode FULL preserves the entire content of the ixpesim output files;
* --stripmode L2 produces a smaller file that is genermane to Level 2 files and
  can be effectively used to test a spectro-polarimetric analysis with the
  photon-list workflow;

Note this requires the input files to be processed with a recent enough gpdsw
version (13.10.0 or later) in order for the conversion to be supported.
"""

from astropy.io import fits
import numpy

from ixpeobssim.core.fitsio import read_hdu_list_in_memory, set_tlbounds
from ixpeobssim.evt.fmt import standard_radec_to_xy, standard_xy_columns_kwargs, \
    set_standard_xy_header_limits, set_object_header_keywords, set_wcs_header_keywords
from ixpeobssim.evt.kislat2015 import xStokesAnalysis
from ixpeobssim.instrument.du import det_name_to_du_id
from ixpeobssim.instrument.gpd import detphi_to_phi
from ixpeobssim.instrument.mma import gpd_to_sky
from ixpeobssim.instrument.sc import pointing_direction
from ixpeobssim.irf.ebounds import channel_to_energy, TLMIN, TLMAX
from ixpeobssim.irfgen.auxiliary import load_pha_model, event_weights, AUX_VERSION
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.math_ import modulo_2pi
from ixpeobssim.utils.os_ import check_output_file

#pylint: disable=invalid-name, too-many-locals, no-member, consider-using-f-string

STRIP_MODES = ('FULL', 'L2')
DEFAULT_STRIP_MODE = 'FULL'


PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_suffix('simfmt')
PARSER.add_weightname()
PARSER.add_auxversion()
PARSER.add_mc()
PARSER.add_argument('--detphiname', type=str, default='DETPHI2',
    help='The column name for the azimuthal angle')
PARSER.add_argument('--stripmode', choices=STRIP_MODES, default=DEFAULT_STRIP_MODE)
PARSER.add_overwrite()


def _remove_wcs_keywords(header, index):
    """Remove all the WCS-related header keywords at a given index.
    """
    logger.debug('Removing WCS header keywords at index %d...', index)
    for key in ('TCTYP', 'TCUNI', 'TCRPX', 'TCRVL', 'TCDLT'):
        key = '%s%d' % (key, index)
        if key not in header:
            logger.warning('Could not find key %s...', key)
        else:
            logger.info('Removing key %s...', key)
            header.remove(key)


def _strip_hdu(hdu, *col_names):
    """Strip down a HDU selecting a subset of the columns.

    Arguments
    ---------
    hdu : FitsHDU instance
        The original HDU to be filtered.

    col_names : iterable
        The column names to be retained in the filtered HDU.
    """
    logger.info('Stripping down %s extension...', hdu.name)
    cols = [col for col in hdu.data.columns if col.name in col_names]
    # At this point we have the additional complication that, depending on
    # whether the SC_DATA extension is present in the original ixpesim file,
    # we might have the WCS information set for the X and Y columns that we just
    # added, and we need to remove them in order not to have *two* WCS defined
    # in the final file. Not even worth mentioning how horrible this hack is.
    header = hdu.header
    if hdu.name == 'EVENTS':
        _col_names = [col.name for col in hdu.data.columns]
        for _col in ('X', 'Y'):
            if _col in _col_names:
                _remove_wcs_keywords(header, _col_names.index(_col) + 1)
    return fits.BinTableHDU.from_columns(cols, header=header)


def _strip_hdu_list_base(hdu_list, config_dict):
    """Strip a list of HDUs based on a configuration dictionary.

    Arguments
    ---------
    hdu_list : FitsHDUList instance
        The original HDU list to be filtered.

    config_dict : dict
        The configuration dictionary, in the form {extension_name: colum_list}
    """
    for ext_name, col_names in config_dict.items():
        hdu_list[ext_name] = _strip_hdu(hdu_list[ext_name], *col_names)


def _strip_hdu_list_l2(hdu_list):
    """Filter an HDU list in the L2 mode.

    Here we only keep the columns of the EVENTS extension that are present in
    Level-2 files and, in addition, we preserve the header of the MONTE_CARLO
    extension, which contains some useful meta-data.
    """
    config_dict = {
        'EVENTS': ('TIME', 'PI', 'W_MOM', 'X', 'Y', 'Q', 'U'),
        'MONTE_CARLO': ()
        }
    _strip_hdu_list_base(hdu_list, config_dict)
    # If the file has been produce through the photon-list mechanism, we can
    # safely remove the SC_DATA extension, which is not necessary in the L2 output.
    for ext_name in ('SC_DATA',):
        try:
            hdu_list.pop(ext_name)
            logger.info('%s extension wiped out...', ext_name)
        except KeyError:
            pass


def format_file(file_path, **kwargs):
    """Format a single file.
    """
    output_file_path = check_output_file(file_path, kwargs.get('suffix'), kwargs.get('overwrite'))
    if output_file_path is None:
        return output_file_path

    # Retrieve the relevant input data.
    hdu_list = read_hdu_list_in_memory(file_path)
    primary_header = hdu_list['PRIMARY'].header
    det_name = primary_header['DETNAM']
    det_id = primary_header['DET_ID']
    try:
        du_id = det_name_to_du_id(det_name)
    except RuntimeError:
        du_id = 1
        logger.warning('DET_ID header keyword is "%s", DU ID set to %d', det_id, du_id)
    # EVENTS extension...
    evt_header = hdu_list['EVENTS'].header
    evt_data = hdu_list['EVENTS'].data
    pha = evt_data['PHA']
    detphi = evt_data[kwargs.get('detphiname')]
    met = evt_data['TIME']
    weights = event_weights(evt_data, kwargs.get('weightname'))
    if kwargs.get('mc'):
        logger.info('Using Monte Carlo absorption point...')
        try:
            mc_data = hdu_list['MONTE_CARLO'].data
        except:
            abort('Input file has no MONTE_CARLO extension, cannot run with --mc True')
        detx = mc_data['ABS_X']
        dety = mc_data['ABS_Y']
    else:
        logger.info('Using reconstructed absorption point...')
        detx = evt_data['ABSX']
        dety = evt_data['ABSY']
    # SC_DATA extension: for files created with ixpesim via the photon list
    # mechanism we should have all the pointing history available in the SC_DATA
    # extension, but files created with the "standard" ixpesim mechanism (e.g.,
    # those created for the purpose of generating response files) will miss this.
    # In this case we set the sc_data local variable to None and we rely on this
    # in the following to detect whether we need to set specific columns or not.
    try:
        sc_header = hdu_list['SC_DATA'].header
        ra_pnt = sc_header['RA_PNT']
        dec_pnt = sc_header['DEC_PNT']
        roll_angle = sc_header['ROLL']
        sc_data = hdu_list['SC_DATA'].data
    except KeyError as e:
        logger.warning(e)
        sc_data = None
        roll_angle = 0.
        logger.warning('Roll angle set to %.3f.', roll_angle)

    # Calculate the missing columns and add them to the EVENTS extension.
    pha_model = load_pha_model(kwargs.get('auxversion', AUX_VERSION))
    pi = pha_model(pha)
    rec_energy = channel_to_energy(pi)
    phi = detphi_to_phi(detphi, du_id, roll_angle)
    # Need to rotate the PHI angle by -90 degrees in order to have the origin
    # of the coordinate system for the position angle at the celestial North,
    # see https://bitbucket.org/ixpesw/ixpeobssim/issues/597
    # There has been a lot of back and forth on this one, and we finally
    # convinced ourselves that we need to rotate by 90 degrees in evt.ixpesim,
    # *and* rotate back by -90 degree here for the whole thing to round-trip
    # correctly, i.e., produce the right pattern of Stokes crosstalk, and
    # preserve the instrinsic source pattern in detector coordinates.
    phi = modulo_2pi(phi - 0.5 * numpy.pi)
    detq = xStokesAnalysis.stokes_q(phi, weights=None)
    detu = xStokesAnalysis.stokes_u(phi, weights=None)
    if sc_data is not None:
        ra, dec = gpd_to_sky(detx, dety, met, *pointing_direction(sc_data, met), du_id, roll_angle)
        x, y = standard_radec_to_xy(ra, dec, ra_pnt, dec_pnt)
        xkwargs, ykwargs = standard_xy_columns_kwargs(ra_pnt, dec_pnt)

    logger.info('Adding necessary columns...')
    cols = hdu_list['EVENTS'].data.columns
    cols += fits.Column(name='DETX', format='E', array=detx)
    cols += fits.Column(name='DETY', format='E', array=dety)
    cols += fits.Column(name='PI', format='E', array=pi)
    cols += fits.Column(name='ENERGY', format='E', array=rec_energy)
    cols += fits.Column(name='PHI', format='E', array=phi)
    if sc_data is not None:
        cols += fits.Column(name='RA', format='E', array=ra)
        cols += fits.Column(name='DEC', format='E', array=dec)
        cols += fits.Column(name='X', format='E', array=x, **xkwargs)
        cols += fits.Column(name='Y', format='E', array=y, **ykwargs)
    cols += fits.Column(name='Q', format='E', array=detq)
    cols += fits.Column(name='U', format='E', array=detu)
    cols += fits.Column(name='W_MOM', format='E', array=weights)
    hdu = fits.BinTableHDU.from_columns(cols, header=evt_header)
    # Now that we have added all the columns, we need to update the EVENTS extension
    # in the underlying HDU list, so that the processing can continue in the
    # proper fashion.
    hdu_list['EVENTS'] = hdu

    # Add the pointing information to the primary header, so that down the road
    # xpbin can make good use of it.
    logger.info('Updating headers...')
    if sc_data is not None:
        set_object_header_keywords(hdu_list['PRIMARY'], ra_pnt, dec_pnt)
    for ext in ('PRIMARY', 'EVENTS', 'GTI'):
        header = hdu_list[ext].header
        header.set('DETNAM', det_name)
        header.set('DET_ID', det_id)

    # Strip the un-necessary columns, if needed.
    if kwargs.get('stripmode') == 'L2':
        _strip_hdu_list_l2(hdu_list)

    # Add TLMIN and TLMAX for the X and Y columns, which is necessary, e.g.,
    # for the file to be properly displayed in ds9. Note this needs to
    # pick up the proper column numbers, and so needs to be done *after* all the
    # relevant columns have been added.
    hdu = hdu_list['EVENTS']
    set_standard_xy_header_limits(hdu)
    set_wcs_header_keywords(hdu)
    # And also, after issue https://github.com/lucabaldini/ixpeobssim/issues/688
    # we need to add the TLMIN and TLMAX header keywords for the PI column.
    set_tlbounds(hdu, 'PI', TLMIN, TLMAX)
    hdu_list['EVENTS'] = hdu

    # Write the processed output file.
    logger.info('Writing processed file to %s...', output_file_path)
    hdu_list.writeto(output_file_path, overwrite=True)
    return output_file_path


def xpsimfmt(**kwargs):
    """Format ground reconstructed data in such a way they are interoperable
    with the ixpeobssim tools---particularly xpbin.
    """
    return [format_file(file_path, **kwargs) for file_path in kwargs.get('filelist')]


def main():
    """main() entry point.
    """
    xpsimfmt(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()

import glob
import os

from astropy.io import fits
from astropy.table import Table
import numpy

from ixpeobssim.core import pipeline
from ixpeobssim.binning.misc import xBinnedLightCurve
from ixpeobssim.evt.gti import xGTIList
from ixpeobssim.evt.fmt import xBinTableHDUGTI
from ixpeobssim.evt.event import xEventFile, xEventFileFriend
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger


PARSER = xArgumentParser()
PARSER.add_argument('folder', type=str,
                    help='path to the input observation folder')

PARSER.add_argument('--l2files', type=str, default=None, nargs='+',
                    help='level 2 file list')


def create_ineclipse_gtis(*hk_file_paths, complement=False):
    """ Create arrays of TSTART and TSTOP corresponding to time periods when 
    the spacecraft is not INSUN (a.k.a. not illuminated by the sun).
    If complement is True, return the INSUN intervals instead"""
    tstart = []
    tstop = []
    for file_path in hk_file_paths:
        hdul = fits.open(file_path)
        data = hdul['HK'].data
        adsec2ecl = data['ADSEC2ECL']
        adsec2sun = data['ADSEC2SUN']
        time = data['TIME']
        if complement:
            mask = (adsec2ecl >= 0) * (adsec2sun < 0)
        else:
            mask = (adsec2ecl < 0) | (adsec2sun >= 0)
        mask = mask.astype(int)
        mask_switch = numpy.ediff1d(mask)
        start_idx = numpy.where(mask_switch > 0)[0] + 1
        if mask[0]:
            start_idx = numpy.append(numpy.array([0]), start_idx)
        stop_idx = numpy.where(mask_switch < 0)[0]
        if mask[-1]:
            stop_idx = numpy.append(stop_idx, numpy.array([len(mask)-1]))
        tstart.append(time[start_idx])
        tstop.append(time[stop_idx])
    return numpy.concatenate(tstart), numpy.concatenate(tstop)
    

def create_gti_extension(start_met, stop_met, tstarts, tstops):
    """ Create a 'GTI' extension from the given arrays of tstarts and tstops.
    This is essentially converting the two arrays of start and stop to 
    an array of tuples (start, stop), as required by the constructor of
    xGTIList."""
    gtis = []
    for start, stop in zip(tstarts, tstops):
        gtis.append([start, stop])
    gti_list = xGTIList(start_met, stop_met, *gtis)
    return xBinTableHDUGTI([gti_list.start_mets(), gti_list.stop_mets()]) 
 

def intersect_obs_gti(obs_file_path, starts, stops):
    """ Create a new 'GTI' extension (of type xBinTableHDUGTI) intersecting the 
    GTIs of the observation with the given starts and stops."""
    obs_file = xEventFile(obs_file_path)
    gti_list = obs_file.get_gti_list()
    new_gtis = gti_list.update(starts, stops)
    return create_gti_extension(obs_file.start_met(), obs_file.stop_met(),
                                *new_gtis)

def create_ineclipse_gti_extension(obs_file_path, *hk_file_paths):
    """ Create a new 'GTI' extension (of type xBinTableHDUGTI) intersecting the 
    GTIs of the observation with the times where the spacecraft is INECLIPSE,
    i.e. not illuminated by the sun.
    """
    inecl_starts, inecl_stops = create_ineclipse_gtis(*hk_file_paths)
    return intersect_obs_gti(obs_file_path, inecl_starts, inecl_stops)


def create_insun_gti_extension(obs_file_path, *hk_file_paths):
    """ Create a new 'GTI' extension (of type xBinTableHDUGTI) intersecting the 
    GTIs of the observation with the times where the spacecraft is INSUN,
    i.e. illuminated by the sun.
    """
    insun_starts, insun_stops = create_ineclipse_gtis(*hk_file_paths,
                                                      complement=True)
    return intersect_obs_gti(obs_file_path, insun_starts, insun_stops)


def gti_mask(obs_file):
    gti_list = obs_file.get_gti_list()
    return gti_list.filter_event_times(obs_file.time_data())


def plot_gtis(gti_data, color='g', alpha=0.3, label=None, **plot_opts):
    """
    """
    starts = gti_data['START']
    stops = gti_data['STOP']
    for start, stop in zip(starts, stops):
        span = plt.axvspan(start, stop, color=color, alpha=alpha, **plot_opts)
    if label is not None:
        span.set_label(label)
        plt.legend()


def write_gti_extension(obs_file_path, gti_extension, tag='inecl'):
    ''' Creates a new file with the GTI table from the inecl/insun GTI 
    extension
    '''
    out_path = os.path.splitext(obs_file_path)[0]+f'_{tag}.fits'
    hdul = fits.open(obs_file_path)
    primary = fits.PrimaryHDU(header = hdul[0].header)
    table1 = fits.BinTableHDU(Table(hdul[1].data))
    table1.header.extend(hdul[1].header, update=True)
    table2 = gti_extension
    new_hdul = fits.HDUList([primary, table1, table2])
    hdul.close()
    new_hdul.writeto(out_path, overwrite=True)
    filter_gtis(out_path)
    return(out_path)


def filter_gtis(file_path):
    ''' Filter a level 2 file according to its GTI table
    '''
    gti = fits.open(file_path)['GTI'].data
    evt = xEventFile(file_path)
    segs = []
    for j in range(len(gti)):
        segs.append(numpy.logical_and((evt.time_data()>gti[j]['START']),
                    (evt.time_data()<gti[j]['STOP'])))
    mask = numpy.logical_or.reduce(segs)
    evt.filter(mask)
    evt.write(file_path, overwrite=True)
    evt.close()


def update_livetime(lvl2_file_path, lvl1_file_path):
    '''
    '''
    friend = xEventFileFriend(lvl2_file_path, lvl1_file_path)
    gti = fits.open(lvl2_file_path)['GTI'].data
    time = friend.l1value('TIME', all_events=True)
    filter = []
    for j in range (len(gti['START'])):
        start = gti['START'][j]
        stop = gti['STOP'][j]
        filter.append(numpy.logical_and((time>start), (time<stop)))
    time_mask = numpy.logical_or.reduce(filter)
    lt_microsec = numpy.sum((friend.l1value('LIVETIME', all_events=True)
                            )*time_mask)
    #gti.close()
    hdul = fits.open(lvl2_file_path, mode='update')
    for hdu in hdul:
        hdu.header['LIVETIME'] = lt_microsec/1.e6
    input (f'new_lt = {lt_microsec/1.e6}')
    hdul.flush()
    hdul.close()


def build_l2_file_dict_from_folder(folder):
    """ Note: this is sorting the files alphabetically, in case this is
    relevant in any way
    """
    l2_folder = os.path.join(folder, 'event_l2')
    logger.debug('Searching for LEVEL 2 files in %s' % l2_folder)
    file_dict = {}
    for du_id in DU_IDS:
        file_list = sorted(glob.glob(
            os.path.join(l2_folder, 'ixpe*_det%s_evt2_v??.fits' % du_id)
        ))
        logger.debug('%d LV2 files found for DU %d: %s' \
                     % (len(file_list), du_id, file_list))
        file_dict[du_id] = file_list
    return file_dict


def build_l2_file_dict_from_file_list(file_list):
    """ Note: this preserve the order of the files given by the user, in case
    this is relevant in any way
    """
    file_dict = {}
    for file_path in file_list:
        if 'det' not in file_path:
            raise ValueError('Cannot find a valid DU number for file %s' % file_path)
        file_name = os.path.basename(file_path)
        pieces = file_name.split('det')
        du_id = int(pieces[1][0])
        if (du_id not in DU_IDS):
            raise ValueError('Cannot find a valid DU number for file %s' % file_path)
        if du_id not in file_dict:
            file_dict[du_id] = [file_path]
        else:
            file_dict[du_id].append(file_path)
    return file_dict


def build_l1_file_dict_from_folder(folder, du_ids=DU_IDS):
    """ Note: this is sorting the files alphabetically, in case this is
    relevant in any way
    """
    l1_folder = os.path.join(folder, 'event_l1')
    logger.debug('Searching for LEVEL 1 files in %s' % l1_folder)
    file_dict = {}
    for du_id in du_ids:
        file_list = sorted(glob.glob(
            os.path.join(l1_folder, 'ixpe*_det%s_evt1_v??.fits' % du_id)
        ))
        logger.debug('%d LV1 files found for DU %d: %s' \
                     % (len(file_list), du_id, file_list))
        file_dict[du_id] = file_list
    return file_dict

def find_hk_file_list_in_folder(folder):
    """
    """
    hk_folder = os.path.join(folder, 'hk')
    logger.debug('Searching for HK files in %s' % hk_folder)
    file_list = sorted(glob.glob(
        os.path.join(hk_folder, 'ixpe*_all_adc_0110_v??.fits')
    ))
    if len(file_list) == 0:
        raise ValueError('Cannot find the appropriate hk file(s) in folder %s' % hk_folder)
    logger.debug('%d hk files found: %s'  % (len(file_list), file_list))
    return file_list


if __name__ == '__main__':
    args = PARSER.parse_args()
    if args.l2files is not None:
        l2_file_dict = build_l2_file_dict_from_file_list(args.l2files)
    else:
        l2_file_dict = build_l2_file_dict_from_folder(args.folder)
    
    l1_file_dict = build_l1_file_dict_from_folder(args.folder,
                                                  du_ids=l2_file_dict.keys())
    hk_file_list = find_hk_file_list_in_folder(args.folder)
    for du_id, l2_file_paths in l2_file_dict.items():
        l1_file_paths = l1_file_dict[du_id]
        for l2_file_path in l2_file_paths:
            logger.info('Processing file %s ...' % l2_file_path)
            ineclipse_gti_ext = create_ineclipse_gti_extension(l2_file_path,
                                                               *hk_file_list)
            inecl_path = write_gti_extension(l2_file_path, ineclipse_gti_ext,
                                             tag='inecl')
            update_livetime(inecl_path, l1_file_paths)
            insun_gti_ext = create_insun_gti_extension(l2_file_path,
                                                       *hk_file_list)
            insun_path = write_gti_extension(l2_file_path, insun_gti_ext,
                                             tag='insun')
            update_livetime(insun_path, l1_file_paths)

            plt.figure()
            plot_gtis(ineclipse_gti_ext.data, label='INECLIPSE')
            plot_gtis(insun_gti_ext.data, color='r', label='INSUN')
            l2_file_path = l2_file_path.replace(' ', '\ ')
            light_curve_path = pipeline.xpbin(l2_file_path, tbins=100000, algorithm='LC', overwrite=True)
            light_curve = xBinnedLightCurve.from_file_list(light_curve_path)
            rate, rate_error = light_curve.rate(), light_curve.rate_error()
            plt.errorbar(light_curve.TIME, rate, rate_error, fmt='o')
            plt.show()

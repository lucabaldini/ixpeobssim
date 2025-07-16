from astropy.io import fits
import matplotlib.patches as mpatches
import numpy as np

from ixpeobssim.core import pipeline
from ixpeobssim.binning.misc import xBinnedLightCurve
from ixpeobssim.evt.gti import xGTIList
from ixpeobssim.evt.fmt import xBinTableHDUGTI
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.matplotlib_ import plt



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
        mask_switch = np.ediff1d(mask)
        start_idx = np.where(mask_switch > 0)[0] + 1
        if mask[0]:
            start_idx = np.append(np.array([0]), start_idx)
        stop_idx = np.where(mask_switch < 0)[0]
        if mask[-1]:
            stop_idx = np.append(stop_idx, np.array([len(mask)-1]))
        tstart.append(time[start_idx])
        tstop.append(time[stop_idx])
    return np.concatenate(tstart), np.concatenate(tstop)
    

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
 

def update_obs_gti(obs_file_path, starts, stops):
    """ Create a new 'GTI' extension (of type xBinTableHDUGTI) intersecting the 
    GTIs of the observation with the given starts and stops.
    If complement is True use the complement of the given starts and stops to
    compute the intersection (see ixpeobssim.evt.gti.xGTIList.complement()). """
    obs_file = xEventFile(obs_file_path)
    gti_extension = obs_file.hdu_list['GTI']
    gti_list = gti_extension.get_gti_list()
    new_gtis = gti_list.update(starts, stops)
    return create_gti_extension(obs_file.start_met(), obs_file.stop_met(),
                                *new_gtis)

def create_ineclipse_gti_extension(obs_file_path, *hk_file_paths):
    """ Create a new 'GTI' extension (of type xBinTableHDUGTI) intersecting the 
    GTIs of the observation with the times where the spacecraft is INECLIPSE,
    i.e. not illuminated by the sun.
    """
    inecl_starts, inecl_stops = create_ineclipse_gtis(*hk_file_paths)
    return update_obs_gti(obs_file_path, inecl_starts, inecl_stops)


def create_insun_gti_extension(obs_file_path, *hk_file_paths):
    """ Create a new 'GTI' extension (of type xBinTableHDUGTI) intersecting the 
    GTIs of the observation with the times where the spacecraft is INSUN,
    i.e. illuminated by the sun.
    """
    insun_starts, insun_stops = create_ineclipse_gtis(*hk_file_paths,
                                                      complement=True)
    return update_obs_gti(obs_file_path, insun_starts, insun_stops)


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


if __name__ == '__main__':
    obs_file_path = '/media/alberto/TOSHIBA EXT/xpe/bkg/01006499/event_l2/ixpe01006499_det2_evt2_v01.fits'
    hk_file_paths = ['/media/alberto/TOSHIBA EXT/xpe/bkg/01006499/hk/ixpe01006401_all_adc_0110_v02.fits',
                    '/media/alberto/TOSHIBA EXT/xpe/bkg/01006499/hk/ixpe01006402_all_adc_0110_v02.fits']
    ineclipse_gti_ext = create_ineclipse_gti_extension(obs_file_path, *hk_file_paths)
    plt.figure()
    plot_gtis(ineclipse_gti_ext.data, label='INECLIPSE')
    insun_gtis = create_insun_gti_extension(obs_file_path, *hk_file_paths)
    plot_gtis(insun_gtis.data, color='r', label='INSUN')
    obs_file_path = obs_file_path.replace(' ', '\ ')
    light_curve_path = pipeline.xpbin(obs_file_path, tbins=10000, algorithm='LC')
    light_curve = xBinnedLightCurve.from_file_list(light_curve_path)
    rate, rate_error = light_curve.rate(), light_curve.rate_error()
    plt.errorbar(light_curve.TIME, rate, rate_error, fmt='o')
    plt.show()


    

from enum import IntEnum
import pathlib

from astropy.io import fits
import matplotlib.patches as mpatches
import numpy as np

from ixpeobssim.core import pipeline
from ixpeobssim.binning.misc import xBinnedLightCurve
from ixpeobssim.evt.gti import xGTIList
from ixpeobssim.evt.fmt import xBinTableHDUGTI
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.matplotlib_ import plt


class EdgeType(IntEnum):
    """
    Enumeration representing the types of edges that can occur in Good Time Intervals (GTIs) 
    and Bad Time Intervals (BTIs).

    Members:
        GTI_START: Beginning of a good time interval.
        GTI_STOP: End of a good time interval.
        BTI_START: Beginning of a bad time interval.
        BTI_STOP: End of a bad time interval.
    """
    GTI_START = 0
    GTI_STOP = 1
    BTI_START = 2
    BTI_STOP = 3


class UnexpectedEdgeType(RuntimeError):
    """
    Custom exception raised when an unexpected sequence of GTI/BTI edges is encountered.

    Attributes:
        msg (str): A detailed error message describing the unexpected edge and current state.
    """
    def __init__(self, edge_time, edge_type, in_old_gti, in_bti):
        """
        Initialize the exception with the edge time, type, and current interval state.

        Args:
            edge_time: The time of the edge event.
            edge_type: The type of the edge (GTI_START, GTI_STOP, etc.).
            in_old_gti (bool): True if currently inside an old GTI.
            in_bti (bool): True if currently inside a BTI.
        """
        self.msg = f'Unexpected edge type {str(EdgeType(edge_type))} at time '\
                   f'{edge_time}, when "in_old_gti" is {in_old_gti} and '\
                   f'"in_bti" is {in_bti}'

    def __str__(self):
        return self.msg


def update_gti(gti_start, gti_stop, bti_start, bti_stop):
    """
    Update Good Time Intervals (GTIs) by removing Bad Time Intervals (BTIs) from them.

    This function merges and processes the GTI and BTI edges in chronological order, producing
    new GTIs that exclude the periods defined as BTIs.

    Args:
        gti_start (np.ndarray): Array of start times for the original GTIs.
        gti_stop (np.ndarray): Array of stop times for the original GTIs.
        bti_start (np.ndarray): Array of start times for BTIs.
        bti_stop (np.ndarray): Array of stop times for BTIs.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Arrays of start and stop times for the new GTIs.
    """
    # Combine all edge times and types
    edge_times = np.hstack((gti_start, gti_stop, bti_start, bti_stop))
    edge_types = np.array(
        [EdgeType.GTI_START] * len(gti_start) +
        [EdgeType.GTI_STOP] * len(gti_stop) +
        [EdgeType.BTI_START] * len(bti_start) +
        [EdgeType.BTI_STOP] * len(bti_stop)
    )

    # Sort edges by time (stable sort ensures start precedes stop if times are equal)
    idx = np.argsort(edge_times, kind='mergesort')
    edge_times = edge_times[idx]
    edge_types = edge_types[idx]

    # State trackers
    in_old_gti = False
    in_bti = False
    new_gti_start = []
    new_gti_stop = []

    # Process each edge in time order
    for edge_time, edge_type in zip(edge_times, edge_types):
        # Case 1: outside both GTI and BTI
        if not in_old_gti and not in_bti:
            if edge_type == EdgeType.GTI_START:
                new_gti_start.append(edge_time)
                in_old_gti = True
            elif edge_type == EdgeType.BTI_START:
                in_bti = True
            else:
                raise UnexpectedEdgeType(edge_time, edge_type, in_old_gti, in_bti)
        # Case 2: inside BTI but not GTI
        elif not in_old_gti and in_bti:
            if edge_type == EdgeType.GTI_START:
                in_old_gti = True
            elif edge_type == EdgeType.BTI_STOP:
                in_bti = False
            else:
                raise UnexpectedEdgeType(edge_time, edge_type, in_old_gti, in_bti)
        # Case 3: inside GTI but not BTI
        elif in_old_gti and not in_bti:
            if edge_type == EdgeType.GTI_STOP:
                in_old_gti = False
                new_gti_stop.append(edge_time)
            elif edge_type == EdgeType.BTI_START:
                in_bti = True
                new_gti_stop.append(edge_time)  # End current GTI at BTI start
            else:
                raise UnexpectedEdgeType(edge_time, edge_type, in_old_gti, in_bti)
        # Case 4: inside both GTI and BTI
        else:
            if edge_type == EdgeType.GTI_STOP:
                in_old_gti = False
            elif edge_type == EdgeType.BTI_STOP:
                in_bti = False
                new_gti_start.append(edge_time)  # Start new GTI after BTI ends
            else:
                raise UnexpectedEdgeType(edge_time, edge_type, in_old_gti, in_bti)

    return np.array(new_gti_start), np.array(new_gti_stop)


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
    gti_data = obs_file.gti_data()
    original_gti_start = gti_data['START']
    original_gti_stop = gti_data['STOP']
    new_gtis = update_gti(original_gti_start, original_gti_stop,
                          starts, stops)
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


    

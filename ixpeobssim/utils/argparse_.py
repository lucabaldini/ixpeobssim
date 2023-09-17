#!/usr/bin/env python
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


import argparse
import ast

from ixpeobssim import IXPEOBSSIM_DATA, IXPEOBSSIM_DOC_FIG_MODELS, IXPEOBSSIM_DOC_FIG_OBSSIM
from ixpeobssim.irf import DEFAULT_IRF_NAME
from ixpeobssim.utils.logging_ import startmsg
from ixpeobssim.utils.time_ import DATETIME_FMT


#pylint: disable = empty-docstring, invalid-name


class xArgumentFormatter(argparse.RawDescriptionHelpFormatter,
                         argparse.ArgumentDefaultsHelpFormatter):

    """Do nothing class combining our favorite formatting for the
    command-line options, i.e., the newlines in the descriptions are
    preserved and, at the same time, the argument defaults are printed
    out when the --help options is passed.

    The inspiration for this is coming from one of the comments in
    https://stackoverflow.com/questions/3853722
    """

    pass



class xArgumentParser(argparse.ArgumentParser):

    """Light-weight wrapper over the argparse ArgumentParser class.

    This is mainly intended to reduce boilerplate code and guarantee a minimum
    uniformity in terms of how the command-line options are expressed across
    different applications.

    Warning
    -------
    Mind you should refrain adding options containing a hyphen, as that will
    break the pipeline code. (This is a consequence of the fact that argparse
    tranforms, for good reasons, hyphens into underscores at parse time.)
    """

    def __init__(self, prog=None, usage=None, description=None):
        """Constructor.
        """
        argparse.ArgumentParser.__init__(self, prog, usage, description,
                                         formatter_class=xArgumentFormatter)

    @staticmethod
    def optional_int(value):
        """Small convenience funtion to convert 'None' into None when parsing
        arguments.
        """
        if value == 'None':
            return None
        return int(value)

    def print_help(self):
        """Overloaded method.

        Here we print the ixpeobssim start message, in addition to the standard
        help output.
        """
        startmsg()
        argparse.ArgumentParser.print_help(self)

    def parse_args(self, args=None, namespace=None):
        """Overloaded method.

        Here we print the ixpeobssim start message, in addition to the standard
        help output.
        """
        args = argparse.ArgumentParser.parse_args(self, args, namespace)
        startmsg()
        return args

    def add_boolean(self, name, default, help):
        """Add a boolean argument.
        """
        self.add_argument(name, type=ast.literal_eval, choices=[True, False],
                          default=default, help=help)

    def add_batch(self):
        """Custom option.
        """
        self.add_argument('--batch', action='store_true', default=False,
                          help='run in batch (do not show plots)')

    def add_configfile(self, default=None, required=True):
        """Custom option.
        """
        if required:
            kwargs = dict(required=True)
        else:
            kwargs = dict(default=None)
        self.add_argument('--configfile', type=str,
                          help='path to the input configuration file',
                          **kwargs)

    def add_deadtime(self, default=0.00108):
        """Custom option.
        """
        self.add_argument('--deadtime', type=float, default=default,
                          help='average dead time per event in s')

    def add_duration(self, default=1000):
        """Custom option.
        """
        self.add_argument('--duration', type=float, default=default,
                          help='duration of the observation in s')

    def add_gti_settings(self, default_min_duration=0., default_start_pad=0., default_stop_pad=0.):
        """Custom option.
        """
        self.add_argument('--gtiminduration', type=float, default=default_min_duration,
                          help='minimum GTI duration in s')
        self.add_argument('--gtistartpad', type=float, default=default_start_pad,
                          help='time padding at the start of each GTI in s')
        self.add_argument('--gtistoppad', type=float, default=default_stop_pad,
                          help='time padding at the stop of each GTI in s')

    def add_ebinning(self, default_emin=2., default_emax=8.,
                     ebin_algs=['FILE', 'LIN', 'LOG', 'LIST'],
                     default_ebinalg='LOG', default_ebins=4):
        """Custom options.
        """
        self.add_ebounds(default_emin, default_emax)
        self.add_argument('--ebinalg', choices=ebin_algs,
                          default=default_ebinalg,
                          help='energy binning specification')
        self.add_argument('--ebins', type=int, default=default_ebins,
                          help='number of bins for LIN/LOG energy binning')
        self.add_argument('--ebinfile', type=str, default=None,
                          help='path to the energy bin definition file')
        self.add_argument('--ebinning', type=ast.literal_eval, default=None,
                          help='list containing the bin edges')

    def add_ebounds(self, default_emin=2., default_emax=8.):
        """Custom options.
        """
        self.add_argument('--emin', type=float, default=default_emin,
                          help='minimum energy in keV')
        self.add_argument('--emax', type=float, default=default_emax,
                          help='maximum energy in keV')

    def add_eef(self, default=1.):
        """Custom option for the encircled energy fraction factor.

        This is a multiplicative factor (less or equal to 1) used in sensitivity
        calculation to account for the effect of the far tails of the PSF, which
        will be removed by any sensible spatial cut to reduce the background.
        """
        self.add_argument('--eef', type=float, default=default,
                          help='encircled-energy-fraction factor')

    def add_file(self):
        """Custom option.
        """
        self.add_argument('file', type=str,
                          help='path to the input file')

    def add_filelist(self):
        """Custom option.
        """
        self.add_argument('filelist', nargs='+',
                          help='path(s) to the input file(s)')

    def add_grayfilter(self):
        """Custom option.
        """
        self.add_boolean('--grayfilter', default=False,
                          help='enable the gray filter')

    def add_irfname(self, default=DEFAULT_IRF_NAME, required=False):
        """Custom option.
        """
        self.add_argument('--irfname', type=str, default=default, required=required,
                          help='name of the response functions to be used')

    def add_mc(self, default=False):
        """Custom option.
        """
        self.add_boolean('--mc', default=default,
                         help='use Monte Carlo information')

    def add_outfile(self, default=None):
        """Custom option.
        """
        self.add_argument('--outfile', type=str, default=default,
                          help='path to the output file')

    def add_outfolder(self, default=IXPEOBSSIM_DATA):
        """Custom option.
        """
        self.add_argument('--outfolder', type=str, default=default,
                          help='path to the output folder')

    def add_overwrite(self, default=True):
        """Custom option.

        Note we are doing this mess with the eval and choices so that
        this can be conveniently wrapped into the pipeline class.
        """
        self.add_boolean('--overwrite', default=default,
                         help='overwrite existing output files')

    def add_phasebounds(self, default_phasemin=None, default_phasemax=None):
        """Custom options.
        """
        self.add_argument('--phasemin', type=float, default=default_phasemin,
                          help='minimum phase for periodic sources')
        self.add_argument('--phasemax', type=float, default=default_phasemax,
                          help='maximum phase for periodic sources')

    def add_roll(self, default=0.):
        """Custom option.
        """
        self.add_argument('--roll', type=float, default=default,
                          help='telescope roll angle in decimal degrees')

    def add_save(self):
        """Custom option.
        """
        self.add_argument('--save', action='store_true', default=False,
                          help='save the output products')

    def add_seed(self, default=None):
        """Custom option.

        The default for the default argument was changed from 0 to None as a
        consequence of issue #189. It is intended that, if the seed is None,
        the downstream code will generate a random random seed instead.
        """
        self.add_argument('--seed', type=self.optional_int, default=default,
                          help='random seed for the simulation')

    def add_srcid(self, default=0):
        """Custom option for the identification of a source into the ROI.
        """
        self.add_argument('--srcid', type=int, default=default,
                          help='the source identifier in the ROITABLE extension')

    def add_srcname(self, default=None):
        """Custom option for the identification of a source into the ROI.
        """
        self.add_argument('--srcname', type=str, default=default,
                          help='the source name in the ROITABLE extension')

    @staticmethod
    def _datefmt():
        """
        """
        fmt_long = DATETIME_FMT.replace('%', '%%')
        fmt_short = fmt_long.split('T')[0]
        return '%s or %s' % (fmt_short, fmt_long)

    def add_startdate(self, default='2022-04-21'):
        """Custom option.

        Note the Easter egg---the default is Martin's birthday!
        """

        self.add_argument('--startdate', type=str, default=default, metavar='DATE',
                          help='observation start date %s' % self._datefmt())

    def add_stopdate(self, default=argparse.SUPPRESS):
        """Custom option.
        """
        self.add_argument('--stopdate', type=str, default=default, metavar='DATE',
                          help='observation stop date %s' % self._datefmt())

    def add_suffix(self, default=None):
        """Custom option.
        """
        self.add_argument('--suffix', type=str, default=default,
                          help='suffix for the output files')

    def add_tbounds(self, default_tmin=None, default_tmax=None):
        """Custom options.
        """
        self.add_argument('--tmin', type=float, default=default_tmin,
                          help='minimum time in s')
        self.add_argument('--tmax', type=float, default=default_tmax,
                          help='maximum tims in s')

    def add_vignetting(self, default=True):
        """
        """
        self.add_boolean('--vignetting', default=default,
                         help='apply MMA vignetting')

    def add_target_source(self):
        """Add the option for the target source.
        """
        group = self.add_mutually_exclusive_group(required=True)
        group.add_argument('--srcname', type=str, default=argparse.SUPPRESS,
                           metavar='NAME',
                           help='name of the target source')
        group.add_argument('--srccoords', type=float, nargs=2,
                           metavar=('RA', 'DEC'), default=argparse.SUPPRESS,
                           help='coordinates of the target source')

    def add_dithering(self, default=True):
        """Add all the dithering-related options.
        """
        self.add_boolean('--dithering', default=default,
                         help='apply the dithering pattern')
        self.add_argument('--ditherampl', type=float, default=1.6,
                          help='the dithering amplitude in arcmin')
        self.add_argument('--ditherpa', type=float, default=907.,
                          help='the dithering main period in s')
        self.add_argument('--ditherpx', type=float, default=101.,
                          help='the dithering x period in s')
        self.add_argument('--ditherpy', type=float, default=449.,
                          help='the dithering y period in s')

    def add_charging(self, default=False):
        """Add all the charging-related options.
        """
        self.add_boolean('--charging', default=default,
                         help='apply the GEM charging')
        self.add_argument('--chrgnside', type=int, default=200,
                          help='number of spatial bins for the charging model')
        self.add_argument('--chrgtstep', type=int, default=30,
                          help='time step for the charging model [s]')
        self.add_argument('--chrgmaps', nargs='+', type=str, default=None,
                          help='maps of the initial charging')
        self.add_argument('--chrgparams', nargs='+', type=str, default=None,
                          help='path to the files of the charging parameters')

    def add_trajectory(self, default=False):
        """Add the trajectory-related options.
        """
        self.add_boolean('--saa', default=default,
                         help='consider the SAA for the GTI calculation')
        self.add_boolean('--occult', default=default,
                         help='consider the Earth occultations for the GTI calculation')

    def add_on_orbit_calibration(self):
        """Add the options related to the onboard calibration.
        """
        self.add_boolean('--onorbitcalib', default=False,
                         help='activate the onboard calibration sources where appropriate')
        self.add_argument('--onorbitcaldemult', type=int, default=1,
                          help='demultiplier for onboard calibration, i.e., do a DU every n orbit(s)')
        self.add_argument('--onorbitcalminduration', type=float, default=600.,
                          help='minimum duration of the onboard calibration runs')
        self.add_argument('--onorbitcalstartpad', type=float, default=30.,
                          help='time padding for starting calibration runs in s')
        self.add_argument('--onorbitcalstoppad', type=float, default=30.,
                          help='time padding for stopping calibration runs in s')
        self.add_argument('--onorbitcalrate', type=float, default=100.,
                          help='average rate in Hz for the onboard Cal C source')

    def add_timeline(self):
        """
        """
        self.add_boolean('--timelinedata', True,
                         help='write the observation timeline into the output file')

    def add_sc_data(self):
        """
        """
        self.add_boolean('--scdata', True,
                         help='write the spacecraft data into the output file')
        self.add_argument('--scdatainterval', type=float, default=5.,
                          help='time interval for the spacecraft data in s')

    def add_phi0(self, default=0.):
        """
        """
        self.add_argument('--phi0', type=float, default=default,
                          help='phase zero value')

    def add_pl_norm(self, default=4.50230e-03):
        """
        """
        self.add_argument('--norm', type=float, default=default,
                          help='the power-law normalization for energy spectrum')

    def add_pl_index(self, default=2.):
        """
        """
        self.add_argument('--index', type=float, default=default,
                          help='the power-law index for energy spectrum')

    def add_column_density(self, default=0.):
        """
        """
        self.add_argument('--nH', type=float, default=default,
                          help='the hydrogen column density (cm^{-2})')

    def add_redshift(self, default=0.):
        """
        """
        self.add_argument('--redshift', type=float, default=default,
                          help='the cosmological redshift')

    def add_stretch(self, default='linear'):
        """
        """
        stretch_values = set(('linear', 'sqrt', 'power', 'log', 'asinh'))
        self.add_argument('--stretch', type=str, default=default,
                          help='The stretch function for image color normalization %s' % stretch_values)

    def add_vrange(self):
        """
        """
        self.add_argument('--vmin', type=float, default=None,
                          help='minimum data range for colormaps')
        self.add_argument('--vmax', type=float, default=None,
                          help='maximum data range for colormaps')

    def add_weights(self, default=False):
        """
        """
        self.add_boolean('--weights', default=default,
                         help='Use weights for the analysis')

    def add_weightname(self, default='alpha075'):
        """
        """
        self.add_argument('--weightname', type=str, default=default,
                          help='name of the weights to be used')

    def add_weightcol(self, default='W_MOM'):
        """
        """
        self.add_argument('--weightcol', type=str, default=default,
                          help='name of the weight column to be used')

    def add_auxversion(self, default=3):
        """
        """
        self.add_argument('--auxversion', type=int, default=default,
                          help='target version for the auxiliary files')

    def add_objname(self, default=None):
        """
        """
        self.add_argument('--objname', type=str, default=None,
                          help='name of the observed object')



class xSourceModelArgumentParser(xArgumentParser):

    """Specialized argument parser for source models.
    """

    def __init__(self, prog=None, usage=None, description=None):
        """Constructor.
        """
        xArgumentParser.__init__(self, prog, usage, description)
        self.add_save()
        self.add_batch()
        self.add_outfolder(default=IXPEOBSSIM_DOC_FIG_MODELS)



class xPipelineParser(xArgumentParser):

    """Specialized argument parser for analysis pipelines.
    """

    def __init__(self, prog=None, usage=None, description=None):
        """Constructor.
        """
        xArgumentParser.__init__(self, prog, usage, description)
        self.add_argument('--target', type=str, default='run',
                          help='the pipeline method to execute')
        self.add_save()
        self.add_batch()
        self.add_outfolder(default=IXPEOBSSIM_DOC_FIG_OBSSIM)
        self.add_overwrite(None)



if __name__ == '__main__':
    parser = xArgumentParser('test', 'toy parser', 'launch me to see options')
    _methods = ['add_argument', 'add_argument_group', 'add_help',
                'add_mutually_exclusive_group', 'add_subparsers']
    for item in dir(parser):
        if item.startswith('add_') and not item in _methods:
            exec('parser.%s()' % item)
    args = parser.parse_args()
    print(args)

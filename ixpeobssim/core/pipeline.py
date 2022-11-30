# Copyright (C) 2016--2022, the ixpeobssim team.
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


"""We need to import the ixpeobssim executables but we don't have an __init__
file in there (and we don't want to add one).

Also, we need to import all the parser object for passing command-line
switches as keyword arguments.
"""

import os
import sys
import glob
import shlex

import numpy

from ixpeobssim import IXPEOBSSIM_BIN, IXPEOBSSIM_DATA, IXPEOBSSIM_DATA_BENCHMARKS
sys.path.append(IXPEOBSSIM_BIN)
sys.dont_write_bytecode = 1

from ixpeobssim.config import config_file_path
from ixpeobssim.binning.polarization import xBinnedPolarizationCube
from ixpeobssim.instrument import DU_IDS, du_suffix
from ixpeobssim.utils.logging_ import logger, abort
import ixpeobssim.utils.matplotlib_ as matplotlib_
import ixpeobssim.utils.os_
from ixpeobssim.utils.os_ import check_input_file

from xpancrkey import xpancrkey as _xpancrkey, PARSER as XPANCRKEY_PARSER
from xpbin import xpbin as _xpbin, PARSER as XPBIN_PARSER
from xpbinview import xpbinview as _xpbinview, PARSER as XPBINVIEW_PARSER
from xpchrgmap import xpchrgmap as _xpchrgmap, PARSER as XPCHRGMAP_PARSER
from xpexposure import xpexposure as _xpexposure, PARSER as XPEXPOSURE_PARSER
from xpgrppha import xpgrppha as _xpgrppha, PARSER as XPGRPPHA_PARSER
from xpmdp import xpmdp as _xpmdp, PARSER as XPMDP_PARSER
from xpobssim import xpobssim as _xpobssim, PARSER as XPOBSSIM_PARSER
from xpophase import xpophase as _xpophase, PARSER as XPOPHASE_PARSER
from xpphase import xpphase as _xpphase, PARSER as XPPHASE_PARSER
from xppicorr import xppicorr as _xppicorr, PARSER as XPPICORR_PARSER
from xpphotonlist import xpphotonlist as _xpphotonlist, PARSER as XPPHOTONLIST_PARSER
from xppimms import xppimms as _xppimms, PARSER as XPPIMMS_PARSER
from xpradialprofile import xpradialprofile as _xpradialprofile, PARSER as XPRADIALPROFILE_PARSER
from xpselect import xpselect as _xpselect, PARSER as XPSELECT_PARSER
from xpsimfmt import xpsimfmt as _xpsimfmt, PARSER as XPSIMFMT_PARSER
from xpsonify import xpsonify as _xpsonify, PARSER as XPSONIFY_PARSER
from xpstripmc import xpstripmc as _xpstripmc, PARSER as XPSTRIPMC_PARSER
from xpstokesalign import xpstokesalign as _xpstokesalign, PARSER as XPSTOKESALIGN_PARSER
from xpstokessmear import xpstokessmear as _xpstokessmear, PARSER as XPSTOKESSMEAR_PARSER
from xpstokesrandom import xpstokesrandom as _xpstokesrandom, PARSER as XPSTOKESRANDOM_PARSER
from xpstokesshuffle import xpstokesshuffle as _xpstokesshuffle, PARSER as XPSTOKESSHUFFLE_PARSER
from xpvisibility import xpvisibility as _xpvisibility, PARSER as XPVISIBILITY_PARSER
from xpxspec import xpxspec as _xpxspec, PARSER as XPXSPEC_PARSER


"""Global setup parameters.

If you're wondering where this is coming from, take a look at
https://en.wikipedia.org/wiki/Run_commands
"""
__rc_params = dict(model='pipeline')


def reset(model_name, **kwargs):
    """
    """
    __rc_params.clear()
    __rc_params.update(model=model_name, **kwargs)


def setup(**kwargs):
    """Setup the global pipeline params.
    """
    __rc_params.update(kwargs)


def set_model(model_name):
    """Set the pipeline model name.

    By default this is the string that will be used, e.g., to create all the
    file path concatenations.
    """
    setup(model=model_name)


def params():
    """Return the underlying rc param dictionary.
    """
    return __rc_params


def param(key, default=None):
    """Retrive a given global configuration parameter.
    """
    return __rc_params.get(key, default)


def model():
    """Return the current model name.
    """
    return param('model')


def save():
    """Return the current value of the 'save' rc parameter.
    """
    return param('save')


def batch():
    """Return the current value of the 'batch' rc parameter.
    """
    return param('batch')


def target():
    """Return the function target for a particular pipeline run.
    """
    return param('target')


def output_folder():
    """Return the current value of the 'outfolder' rc parameter.
    """
    return param('outfolder')


def overwrite():
    """Return the current value of the 'overwrite' rc parameter.
    """
    return param('overwrite')


def figure_name(name):
    """Small convience function to enfore a minimum uniformity in the
    naming scheme for the output file.

    This will prepeng the model name to any name the use passes as an
    argument, and change the latter to all lower case.
    """
    return '%s_%s' % (model(), name.lower())


def set_gcf_name(name):
    """Set the canvas title for the current figure, following the same rules
    of the figure_name() method.
    """
    matplotlib_.plt.gcf().canvas.set_window_title(figure_name(name))


def figure(name):
    """Return a matplotlib figure with the model name prepended to the
    actual name.
    """
    return matplotlib_.plt.figure(figure_name(name))


def residual_figure(name):
    """Create a figure for residual plots.
    """
    return matplotlib_.residual_plot(figure_name(name))


def bootstrap_pipeline(model_name):
    """Convenience bootstrap function.
    """
    from ixpeobssim.utils.argparse_ import xPipelineParser
    parser = xPipelineParser()
    args = parser.parse_args()
    logger.info('Bootstrapping analysis pipeline...')
    reset(model_name, target=args.target)
    setup(**args.__dict__)
    logger.info('Done, %s.' % params())
    import __main__
    try:
        function = __main__.__getattribute__(args.target)
    except AttributeError:
        abort('Cannot execute pipeline %s() target' % args.target)
    function()
    if save():
        matplotlib_.save_all_figures(output_folder(), ['png', 'pdf'])
    if not batch():
        matplotlib_.plt.show()


def suffix(label=None, index=None):
    """Create a suffix to be appended to file names when the pipeline tools
    are created.
    """
    if label is None and index is None:
        return ''
    if index is None:
        return label
    return '%s%04d' % (label, index)


def file_list(*args, **kwargs):
    """Create a file list from a series of patterns.

    This method is very handy to retrieve the files generated by any given
    step of a pipeline in order to process them in the following step, and takes
    care automagically of looping over the three detector units.

    The arguments can be either a string or a (str, int) tuple, in which case
    the string is a label and the integer is a uid attached to it (e.g., the
    identifier for a particular phase selection for a periodic source).
    The arguments are concatenated in the file name with an "_" character.

    Assuming that we have a pipieline bootstrapped with a model called "mymodel",
    the resolution rules will yield:
    pipeline.file_list() -> '$IXPEOBSSIM_DATA/mymodel_du1/2/3.fits'
    pipeline.file_list(('sel', 1)) -> '$IXPEOBSSIM_DATA/mymodel_du1/2/3_sel0001.fits'
    pipeline.file_list(('sel', 1), 'cmap') -> '$IXPEOBSSIM_DATA/mymodel_du1/2/3_sel0001_cmap.fits'

    The label "pha1*" is peculiar in that, if encountered it is automatically
    expanded to include all the flavors of the Stokes parameter spectra.

    Note that, in order to maintain compatibility with Python 2, we opted for
    the current function signature over what would have been the most
    natural way to do things in Python 3, i.e.,
    `file_list(*args, folder_path=IXPEOBSSIM_DATA, check_files=True)`
    """
    # Recursive hack for the "pha1*" label...
    try:
        index = args.index('pha1*')
        _file_list = []
        for alg in ('pha1', 'pha1q', 'pha1u'):
            _args = list(args)
            _args[index] = alg
            _file_list += file_list(*_args, **kwargs)
        return _file_list
    except ValueError:
        pass
    # Main function body.
    folder_path = kwargs.get('folder_path', IXPEOBSSIM_DATA)
    check_files = kwargs.get('check_files', True)
    _file_list = []
    for du_id in DU_IDS:
        file_name = '%s_%s' % (model(), du_suffix(du_id))
        for arg in args:
            if arg is None:
                pass
            elif isinstance(arg, str):
                file_name += '_%s' % suffix(arg)
            else:
                file_name += '_%s' % suffix(*arg)
        file_name = '%s.fits' % file_name
        file_path = os.path.join(folder_path, file_name)
        if check_files:
            check_input_file(file_path, 'fits')
        _file_list.append(file_path)
    return _file_list


def _command_line_switches(*args, **kwargs):
    """Turn a set of positional and/or keyword arguments into a list of
    command-line switches (i.e., strings) that is understood by argparse.

    This will turn, e.g.,

    args = ['file1', 'file2']
    kwargs = {'configfile': 'test.py', 'duration': 100}

    into

    ['file1', 'file2', '--configfile', 'test.py', '--duration', '100']

    Note that the global overwrite flag (defaulting to False) is peculiar in
    that is ruling over the local values passed via command-line arguments to
    the wrapper functions running the ixpeobssim tools.
    """
    # Control the global overwrite pipeline setting.
    if param('overwrite') is not None:
        kwargs.update(overwrite=param('overwrite'))
    if not 'overwrite' in kwargs:
        kwargs.update(overwrite=False)
    switches = ''
    # Loop over the positional arguments first.
    for arg in args:
        switches += ' %s' % arg
    # Add an extra space, if needed, to separate positional arguments from
    # command-line arguments.
    if len(args):
        switches += ' '
    # And now the keyword arguments.
    for key, value in kwargs.items():
        # Need some extra care for lists...
        if isinstance(value, list):
            value = ('%s' % value).replace(' ', '')
        # and tuples...
        if isinstance(value, tuple):
            value = str(value).replace(' ', '')
        # ... and numpy arrays.
        if isinstance(value, numpy.ndarray):
            value = numpy.array2string(value, separator=',')
            value = value.replace(' ', '')
        # Horrible hack for the --modelfiles switch in xpstokesalign---here we
        # are passing two files and we need a way to parse the input in a way
        # that makes the pipeline happy...
        # Note that, in this case, we don't want the file paths to be quoted,
        # so the addition of quotes for strings containing spaces is done in the
        # else branch of this conditional statement.
        # And, for completeness, the same holds for --chrgmaps, --chrgparams,
        # --cmapfiles and --arffiles
        if key in ('modelfiles', 'chrgmaps', 'chrgparams', 'cmapfiles', 'arffiles'):
            if value is not None:
                value = value.strip('[]').replace(',', ' ').replace('\'', '')
        else:
            # And here we try and deal with string command-line options values
            # containing spaces---they need to be quoted.
            if isinstance(value, str) and ' ' in value:
                value = '"%s"' % value
        # Another hack to handle XSPEC models, which are string that possibly
        # contain spaces.
        if key in ['model', 'specmodel', 'polmodel']:
            value = value.replace(' ', '')
        # Another nice one: argparse doesn't seem able to understand
        # negative numbers in engineering notation...
        # Here we need to quote and to add an extra space.
        if isinstance(value, float) and value < 0:
            value = '" %s"' % value
        switches += '--%s %s ' % (key, value)
    # Strip and split. Note that we use the split function from the shlex
    # module, as we want to preserve the quoted arguments, if any, see
    # https://docs.python.org/3/library/shlex.html
    switches.strip()
    return shlex.split(switches)


def _parse_args(parser, *args, **kwargs):
    """Run a set of keyword arguments through a specific argument parser and
    return a full set of options, including all the default values from
    the parser itsels.
    """
    switches = _command_line_switches(*args, **kwargs)
    return parser.parse_args(switches).__dict__


def _update_configfile_kwarg(**kwargs):
    """Small utility function to update a set of command-line option with the
    `configfile` key.

    This is used in all the applications accepting (o requiring) a `configfile`
    command-line switch, so that it can be omitted in a pipeline context where
    the model name is known and the file can be picked up automatically from the
    config folder.
    """
    if not 'configfile' in kwargs:
        kwargs.update(configfile=config_file_path(model()))
    return kwargs


def xpancrkey(*args, **kwargs):
    """App wrapper.
    """
    return _xpancrkey(**_parse_args(XPANCRKEY_PARSER, *args, **kwargs))


def xpbin(*args, **kwargs):
    """App wrapper.
    """
    return _xpbin(**_parse_args(XPBIN_PARSER, *args, **kwargs))


def xpbinview(*args, **kwargs):
    """App wrapper.
    """
    return _xpbinview(**_parse_args(XPBINVIEW_PARSER, *args, **kwargs))


def xpchrgmap(*args, **kwargs):
    """App wrapper.
    """
    return _xpchrgmap(**_parse_args(XPCHRGMAP_PARSER, *args, **kwargs))


def xpexposure(*args, **kwargs):
    """App wrapper.
    """
    return _xpexposure(**_parse_args(XPEXPOSURE_PARSER, *args, **kwargs))


def xpgrppha(*args, **kwargs):
    """App wrapper.
    """
    return _xpgrppha(**_parse_args(XPGRPPHA_PARSER, *args, **kwargs))


def xpmdp(**kwargs):
    """App wrapper.
    """
    kwargs = _update_configfile_kwarg(**kwargs)
    return _xpmdp(**_parse_args(XPMDP_PARSER, **kwargs))


def xpobssim(**kwargs):
    """App wrapper.
    """
    kwargs = _update_configfile_kwarg(**kwargs)
    return _xpobssim(**_parse_args(XPOBSSIM_PARSER, **kwargs))


def xpophase(*args, **kwargs):
    """App wrapper.
    """
    return _xpophase(**_parse_args(XPOPHASE_PARSER, *args, **kwargs))


def xpphase(*args, **kwargs):
    """App wrapper.
    """
    return _xpphase(**_parse_args(XPPHASE_PARSER, *args, **kwargs))


def xppicorr(*args, **kwargs):
    """App wrapper.
    """
    return _xppicorr(**_parse_args(XPPICORR_PARSER, *args, **kwargs))


def xpphotonlist(*args, **kwargs):
    """App wrapper.
    """
    kwargs = _update_configfile_kwarg(**kwargs)
    return _xpphotonlist(**_parse_args(XPPHOTONLIST_PARSER, *args, **kwargs))


def xppimms(**kwargs):
    """App wrapper.
    """
    return _xppimms(**_parse_args(XPPIMMS_PARSER, **kwargs))


def xpradialprofile(**kwargs):
    """App wrapper.
    """
    return _xpradialprofile(**_parse_args(XPRADIALPROFILE_PARSER, **kwargs))


def xpselect(*args, **kwargs):
    """App wrapper.
    """
    return _xpselect(**_parse_args(XPSELECT_PARSER, *args, **kwargs))


def xpsimfmt(*args, **kwargs):
    """App wrapper.
    """
    return _xpsimfmt(**_parse_args(XPSIMFMT_PARSER, *args, **kwargs))


def xpsonify(*args, **kwargs):
    """App wrapper.
    """
    return _xpsonify(**_parse_args(XPSONIFY_PARSER, *args, **kwargs))


def xpstokesalign(*args, **kwargs):
    """App wrapper.
    """
    return _xpstokesalign(**_parse_args(XPSTOKESALIGN_PARSER, *args, **kwargs))


def xpstokessmear(*args, **kwargs):
    """App wrapper.
    """
    return _xpstokessmear(**_parse_args(XPSTOKESSMEAR_PARSER, *args, **kwargs))


def xpstokesrandom(*args, **kwargs):
    """App wrapper.
    """
    return _xpstokesrandom(**_parse_args(XPSTOKESRANDOM_PARSER, *args, **kwargs))


def xpstokesshuffle(*args, **kwargs):
    """App wrapper.
    """
    return _xpstokesshuffle(**_parse_args(XPSTOKESSHUFFLE_PARSER, *args, **kwargs))


def xpstripmc(*args, **kwargs):
    """App wrapper.
    """
    return _xpstripmc(**_parse_args(XPSTRIPMC_PARSER, *args, **kwargs))


def xpvisibility(*args, **kwargs):
    """App wrapper.
    """
    return _xpvisibility(**_parse_args(XPVISIBILITY_PARSER, *args, **kwargs))


def xpxspec(*args, **kwargs):
    """App wrapper.
    """
    return _xpxspec(**_parse_args(XPXSPEC_PARSER, *args, **kwargs))


def standard_ensamble_processing(file_list, mc=True):
    """Standard processing routine.
    """
    # Generate the Stokes spectra.
    for algorithm in ['PHA1', 'PHA1Q', 'PHA1U', 'PHA1QN', 'PHA1UN']:
        xpbin(*file_list, algorithm=algorithm)
    # Generate the modulation cubes, in both the standard and MC flavors.
    kwargs = dict(emin=2., emax=8., ebins=1)
    xpbin(*file_list, algorithm='PCUBE', **kwargs)
    if mc:
        xpbin(*file_list, algorithm='PCUBE', mc=True, suffix='mc_pcube', **kwargs)


def generate_ensamble(size=10, start_seed=0, processing_function=standard_ensamble_processing,
                      cleanup=True, **kwargs):
    """Generate an ensable of simulations to study the a-posteriory statistical
    properties of the measurement setup.

    Arguments
    ---------
    size : int
        The number of independent realizations to be generated

    start_seed : int
        The first random seed to be used (will be incremented by one at each step)

    cleanup : bool
        If True, remove the (potentially large) event files after the processing.

    **Kwargs:
        All the keyword arguments to be passes to xpobssim.
    """
    for key in ['outfile', 'configfile', 'seed']:
        if key in kwargs:
            logger.error('Argument "%s" to pipeline.generate_ensamble() will be ignored!', key)
            kwargs.pop(key)
    for seed in range(start_seed, start_seed + size):
        # Assemble the path to the output file.
        file_name = '%s_rnd%04d' % (model(), seed)
        file_path = os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, file_name)
        # Run the simulation.
        file_list = xpobssim(seed=seed, outfile=file_path, **kwargs)
        # Call the single-realization processing routine.
        if processing_function is not None:
            processing_function(file_list)
        # Remove the (large) event files.
        if cleanup:
            for file_path in file_list:
                ixpeobssim.utils.os_.rm(file_path)


def glob_ensamble(seed, *patterns):
    """Glob an ensamble folder searching for the files matching a particular
    seed and pattern.

    Note this is completely different in spirit wrt file_list(), as the list
    is not assembled a priori, but globbed directly from the file system (i.e.,
    you are guaranteed that the files in the list do exist).
    """
    file_list = []
    for pattern in patterns:
        file_name = '%s_rnd%04d_du?_%s.fits' % (model(), seed, pattern)
        file_list += glob.glob(os.path.join(IXPEOBSSIM_DATA_BENCHMARKS, file_name))
    if len(file_list) % 3 != 0:
        msg = 'glob_ensamble() returning %d file(s)---not a multiple of 3.'
        logger.error(msg, len(file_list))
    # Nice to have the thing sorted...
    file_list.sort()
    return file_list


def glob_ensamble_stokes_spectra(seed, normalized=False):
    """Specialized function to get hold of a consistent set of Stokes spectra
    for a given seed.

    If the normalized argument is False, this returns the PHA1, PHA1Q and PHA1U
    files for the three detector units, otherwise the PHA1QN and PHA1UN spectra
    are returned.
    """
    if not normalized:
        return glob_ensamble(seed, 'pha1', 'pha1q', 'pha1u')
    return glob_ensamble(seed, 'pha1qn', 'pha1un')


def fit_ensamble_stokes_spectra(seed, normalized=False, **kwargs):
    """Fit an ensamble file with XSPEC.

    This is a small convenience function that, under the hood, collects the
    relevant file list, runs xpxspec, and return a list with the best-fit
    parameters and errors.
    """
    kwargs.setdefault('plot', False)
    file_list = glob_ensamble_stokes_spectra(seed, normalized)
    fit_output = xpxspec(*file_list, **kwargs)
    values = []
    for _, value, error, _, _, _ in fit_output:
        values += [value, error]
    return values


def post_process_ensamble_pcubes(seed, label='pcube'):
    """Post-process the polarization cubes for a given ensamble run.
    """
    file_list = glob_ensamble(seed, label)
    cube = xBinnedPolarizationCube.from_file_list(file_list)
    return cube.polarization()

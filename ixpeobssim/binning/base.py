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

"""Base classes for binned data products.
"""

from __future__ import print_function, division

from collections.abc import Iterable
from functools import reduce

import numpy

from ixpeobssim.core.fitsio import xPrimaryHDU, read_hdu_list_in_memory
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.irf import load_arf
from ixpeobssim.utils.astro import build_wcs
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.units_ import arcmin_to_arcsec, arcsec_to_degrees

# pylint: disable=invalid-name, too-many-arguments, no-member


class xEventBinningBase:

    """Base class for the event binning.

    This is essentially opening an event file and keeping track of the
    xEventFile object, along all the keyword arguments passed to the
    constructor.
    """

    INTENT = None
    SUPPORTED_KWARGS = None

    def __init__(self, file_path, **kwargs):
        """Constructor.
        """
        self.event_file = xEventFile(file_path)
        self.irf_name = kwargs.get('irfname')
        if self.irf_name is None:
            self.irf_name = self.event_file.irf_name()
        self.kwargs = kwargs
        self.process_kwargs()

    def process_kwargs(self):
        """Check the keyword arguments and if the output file is not set, create
        a path on the fly, based on the event file and the binning algorithm.
        """
        if self.get('suffix') is None:
            suffix = self.get('algorithm').lower()
        else:
            suffix = self.get('suffix')
        evfile = self.event_file.file_path()
        outfile = evfile.replace('.fits', '_%s.fits' % suffix)
        self.set('outfile', outfile)

    def process_image_ref_kwargs(self):
        """Set the xref and yref values for the output image (in sky) coordinates.

        If either xref or yref are not specifified, the center of the ROI specified
        in the proper extension is assumed.
        """
        # Retieve the coordinates of the ROI center.
        xref, yref = self.event_file.wcs_reference()
        if self.get('xref') is None:
            self.set('xref', xref)
        if self.get('yref') is None:
            self.set('yref', yref)
        # And if we still haven't succeded, we give up.
        if None in (self.get('xref'), self.get('yref')):
            logger.error('Cannot retrieve the ROI center information.')
            abort('Please set the --xref and --yref command-line switches explicitely')

    @staticmethod
    def _image_wcs_kwargs():
        """Return the list of valid keywords for the _build_image_wcs() method.
        """
        return ['npix', 'pixsize', 'proj', 'xref', 'yref']

    @staticmethod
    def _build_image_wcs(default_img_side=17.5, **kwargs):
        """Build the WCS object for the output image.

        This is a convenicence function meant to be used by all the binning
        algorithms needing to create programmatically a WCS object to operate
        in sky coordinates.

        Arguments
        ---------
        default_img_side : float
            The default size of the output image in arcmin.

        kwargs : dict
            All the command-line keyword arguments passed to xpbin.
        """
        ra0 = kwargs.get('xref', 0.)
        dec0 = kwargs.get('yref', 0.)
        nside = kwargs.get('npix', 200)
        # Retrieve the pixel size...
        pixel_size = kwargs.get('pixsize', None)
        # ...fall back to the default image size if the pixel size is not set.
        # (Mind the image size in arcmin, and the pixel size is in arcsec.)
        if pixel_size is None:
            pixel_size = arcmin_to_arcsec(default_img_side) / nside
        # Convert from arcsec to degrees.
        pixel_size = arcsec_to_degrees(pixel_size)
        projection = kwargs.get('proj', 'TAN')
        return build_wcs(ra0, dec0, nside, pixel_size, projection)

    @staticmethod
    def _pixelize_skycoords(ra, dec, wcs_, origin=0, offset=0.5):
        """Convert sky-coordinates into pixel coordinates for a given wcs.

        This is supposed to be the main function for handling the sky-to-pixel
        conversion in binned object, and is fundamentally different from the
        astropy.WCS.wcs_world2pix() method called under the hood, in that
        it provides a series of digitized values suitable to be passed to the
        numpy.histogramdd() machinery, along with the proper binning to be
        used downstream. For this very purpose, an offset in the (0, 1) open
        interval (0.5 by default) is added to the return values, so that the
        x and y output values do not lie on the bin edges when used to create
        histograms downstream. Note that if the offset is set to zero, this will
        introduce an offset in the output maps, see the test_wcs.py unit test
        for more information.

        Note that the origin argument is largely irrelevant as long as the proper
        binning (i.e., the one returned by the function) is used downstream, see
        https://docs.astropy.org/en/stable/api/astropy.wcs.WCS.html#astropy.wcs.WCS.wcs_world2pix
        for more information about the origin.

        Arguments
        ---------
        ra : array_like
            The array of input ra values.

        dec : array_like
            The array of input dec values.

        wcs_ : astropy.wcs.WCS opject
            The WCS object handling the conversion behind the scenes.

        origin : 0 or 1
            The origin for the conversion. See
            https://docs.astropy.org/en/stable/api/astropy.wcs.WCS.html#astropy.wcs.WCS.wcs_world2pix
            for more details.

        offset : float
            Constant offset to be added to the output x and y arrays.
        """
        # Pack the input coordinates in the proper form for the calculation.
        coords = numpy.vstack((ra, dec)).transpose()
        # Call the underlying astropy.WCS.wcs_world2pix() method.
        pix = wcs_.wcs_world2pix(coords, origin)
        # Unpack the output of the conversion---note that we are swapping the
        # axis to match the FITS convension, since this is going to be stored
        # into a FITS file. This avoid a transposition at the read-back time.
        x = pix[:, 1] + offset
        y = pix[:, 0] + offset
        # Since we have the wcs *and* the origin used to do the world-to-pixel
        # conversion, this is *the* right place to construct the binning to
        # be used along with x and y downstream.
        nx, ny = wcs_.array_shape
        xbinning = numpy.linspace(origin, nx + origin, nx + 1)
        ybinning = numpy.linspace(origin, ny + origin, ny + 1)
        return x, y, (xbinning, ybinning)

    def get(self, key, default=None):
        """Convenience method to address the keyword aguments.
        """
        return self.kwargs.get(key, default)

    def set(self, key, value):
        """Convenience method to set keyword arguments.
        """
        logger.info('Setting %s to %s...', key, value)
        self.kwargs[key] = value

    def weight_data(self):
        """Retrieve the weights from the underlying event file.

        This encapsulates the simple logic used downstream by the binning algorithms
        supporting weights, i.e., if the command-line argument `weights` is False
        the weights are all one, while if it is True returns they are the content
        of the column indicated by the `weightcol` command-line argument.
        """
        if not self.get('weights'):
            logger.info('Performing un-weighted analysis...')
            return numpy.full(self.event_file.num_events(), 1.)
        col_name = self.get('weightcol')
        logger.info('Retrieving weights from column "%s"...', col_name)
        return self.event_file.event_data[col_name]

    def check_pcube_weighting_scheme(self, aeff):
        """Simple check on the weighting scheme being used.

        For more background information, see
        https://bitbucket.org/ixpesw/ixpeobssim/issues/573 and
        https://bitbucket.org/ixpesw/ixpeobssim/issues/613
        """
        weight_scheme = aeff.weighting_scheme()
        if self.get('weights') and weight_scheme != 'SIMPLE':
            logger.error('Wrong arf weighting scheme (%s) for %s with options %s',
                weight_scheme, self.__class__.__name__, self.kwargs)
            logger.info('See https://bitbucket.org/ixpesw/ixpeobssim/issues/573 for more info.')
            raise RuntimeError('Cannot create binned data product')

    def load_aeff_for_polarization_analysis(self):
        """Load the proper arf file for a model-independent polarization analysis
        (i.e., to be used for polarization cubes, maps and map cubes).

        This is loading the proper arf file making sure that, when weights are
        used, the SIMPLE weighting prescription is picked.
        """
        aeff = load_arf(self.irf_name, self.event_file.du_id(),
            simple_weighting=self.get('weights'), gray_filter=self.get('grayfilter'))
        self.check_pcube_weighting_scheme(aeff)
        return aeff

    def build_primary_hdu(self, data=None, header=None):
        """Build the primary HDU for the output file.
        """
        primary_hdu = xPrimaryHDU(data, header, 'xpbin.py')
        primary_hdu.setup_header(self.event_file.primary_keywords())
        primary_hdu.add_keyword('BINALG', self.get('algorithm'),
                                'the binning algorithm used')
        primary_hdu.add_comment('%s run with kwargs %s' %\
                                (self.__class__.__name__, self.kwargs))
        return primary_hdu

    @staticmethod
    def read_binning_from_file(file_path):
        """Read a custom binning from file and return a numpy array.
        """
        return numpy.loadtxt(file_path)

    @staticmethod
    def equipopulated_binning(num_bins, data, min_val=None, max_val=None):
        """Create an equipopulated binning based on the values of a generic
        data column.

        Arguments
        ---------
        num_bins : int
            The number of bins

        data : array
            The underlying data to be used for the binning

        min_val : float (optional)
            The minimum data value to be used for the binning

        max_val : float (optional)
            The maximum data value to be used for the binning
        """
        if min_val is None:
            min_val = data.min()
        if max_val is None:
            max_val = data.max()
        mask = (data >= min_val) * (data <= max_val)
        data = data[mask]
        data.sort()
        binning = [min_val]
        for i in range(1, num_bins):
            index = int(i * len(data) / float(num_bins))
            binning.append(data[index])
        binning.append(max_val)
        return numpy.array(binning)

    # pylint: disable=inconsistent-return-statements
    @staticmethod
    def make_binning(bin_alg, min_val=None, max_val=None, num_bins=None,
        bin_list=None, bin_data=None, bin_file=None):
        """Generic function to define a binning with several possible
        different algorithms.
        """
        assert bin_alg in ['LIN', 'LOG', 'LIST', 'EQP', 'FILE']
        if bin_alg == 'LIN':
            assert min_val is not None
            assert max_val is not None
            assert num_bins is not None
            return numpy.linspace(min_val, max_val, num_bins + 1)
        if bin_alg == 'LOG':
            assert min_val is not None and min_val > 0.
            assert max_val is not None
            assert num_bins is not None
            return numpy.logspace(numpy.log10(min_val), numpy.log10(max_val), num_bins + 1)
        if bin_alg == 'EQP':
            assert num_bins is not None
            assert bin_data is not None
            return xEventBinningBase.equipopulated_binning(num_bins, bin_data, min_val, max_val)
        if bin_alg == 'FILE':
            assert bin_file is not None
            return xEventBinningBase.read_binning_from_file(bin_file)
        if bin_alg == 'LIST':
            assert bin_list is not None
            return numpy.array(bin_list, 'd')

    @staticmethod
    def _energy_binning_kwargs():
        """Return the list of valid keywords for the _build_image_wcs() method.
        """
        return ['emin', 'emax', 'ebins', 'ebinalg', 'ebinning', 'ebinfile']

    @staticmethod
    def make_energy_binning(bin_data=None, **kwargs):
        """Small convenience function for energy binning.
        """
        args = (kwargs.get('ebinalg'), kwargs.get('emin'), kwargs.get('emax'),
                kwargs.get('ebins'), kwargs.get('ebinning'), bin_data,
                kwargs.get('ebinfile'))
        binning = xEventBinningBase.make_binning(*args)
        logger.info('Energy binning: %s', binning)
        return binning

    @staticmethod
    def bin_centers(bin_edges):
        """Return an array of bin centers given an array of bin edges.

        Arguments
        ---------
        bin_edges : 1-d array of length (n + 1).
            The array with the bin edges.

        Returns
        -------
        1-d array of length n.
            The array with the values of the bin centers.
        """
        assert bin_edges.ndim == 1
        return 0.5 * (bin_edges[:-1] + bin_edges[1:])

    @staticmethod
    def bin_widths(bin_edges):
        """Return an array of bin widths given an array of bin edges.

        Arguments
        ---------
        bin_edges : 1-d array of length (n + 1).
            The array with the bin edges.

        Returns
        -------
        1-d array of length n.
            The array with the values of the bin widths.
        """
        assert bin_edges.ndim == 1
        return bin_edges[1:] - bin_edges[:-1]

    def bin_(self):
        """Do-nothing method to be reimplemented in the derived classes.
        """
        raise NotImplementedError

    def write_output_file(self, hdu_list, overwrite=True):
        """Write the binned data to the output file.
        """
        hdu_list.info()
        logger.info('Writing %s binned data to %s...', self.get('algorithm'),
                    self.get('outfile'))
        hdu_list.writeto(self.get('outfile'), overwrite=overwrite)
        logger.info('Done.')



class xBinnedFileBase:

    """Base class for binned files.
    """

    # pylint: disable=protected-access

    def __init__(self, file_path):
        """Constructor.

        This is the basic interface to binned FITS files written by xpbin.
        There are two different ways to instantiate class object: (i) call the
        constructor given a file path; (ii) call the class method
        from_file_list() passing a list of file paths (in which case the
        contents of the different files are summed up together according to the
        specific implementation of the __add__() method for any particular
        class).

        Derived classes are responsible for reimplementing the _read_data()
        method, where all the relevant information is extracted from the FITS
        file itself. This was originally implemented because the underlying
        FITS file was closed right away in the constructor, but this prevented
        any easy implementation of a method to save to disk any binned data
        product. So the FITS file is now open for the entire lifetime of the
        object, but the machinery for the internal storage of the data is
        retained, as is now pervasive in the code and it would be too disruptive
        to change it.
        """
        # Mind we Initialize this to avoid isses in the destructor, where we
        # clean things up.
        self.__data_dict = {}
        self.hdu_list = read_hdu_list_in_memory(file_path)
        self.hdu_list.info()
        self.primary_header = self.hdu_list['PRIMARY'].header
        self._read_data()

    @classmethod
    def from_file_list(cls, file_list):
        """Method to instance the binned class from a list of files.
        """
        assert isinstance(file_list, Iterable)
        return reduce(cls.__iadd__, [cls(file_path) for file_path in file_list])

    def _read_data(self):
        """Do-nothing method to be reimplemented in the derived classes.
        """

    def _read_binary_table_data(self, extension_name):
        """Convenience function retrieving all the data from a binary table.

        This is looping over all the columns of the binary table identified
        by a given extension name and inderting all the corresponding arrays
        into the self.__data_dict dictionary.
        """
        data = self.hdu_list[extension_name].data
        for col in data.columns:
            self.__data_dict[col.name] = col.array

    def _read_image_data(self, extension_name):
        """Convenience function retrieving data from an image.
        """
        self.__data_dict[extension_name] = self.hdu_list[extension_name].data

    def _weighted_average(self, other, col_name, weights_col_name, default=0., invert_w2=False):
        """Convenience function to calculate the weighted average of two
        different binned object for a particular column.

        This comes handy when summing binned files, e.g., for the average
        energy or effective modulation factor.

        Note we make all efforts to avoid any zero-division error.
        """
        v1 = self.__data_dict[col_name]
        w1 = self.__data_dict[weights_col_name]
        mask1 = w1 > 0.
        v2 = other.__data_dict[col_name]
        w2 = other.__data_dict[weights_col_name]
        mask2 = w2 > 0.
        if invert_w2:
            w2 = -w2
        # Initialize the array holding the weighted average to the default value.
        val = numpy.full(v1.shape, default)
        # Both weights are non-zero---we can do the actual weighted average.
        mask = numpy.logical_and(mask1, mask2)
        val[mask] = (v1 * w1 + v2 * w2)[mask] / (w1 + w2)[mask]
        # Weights for other are zero.
        mask = numpy.logical_and(mask1, numpy.logical_not(mask2))
        val[mask] = v1[mask]
        # Weights for self are zero.
        mask = numpy.logical_and(numpy.logical_not(mask1), mask2)
        val[mask] = v2[mask]
        assert numpy.isfinite(val).all()
        return val

    def __getattr__(self, name):
        """Proxy method to retrieve the underlying data as if they were class
        members.

        Warning
        -------
        I guess I was feeling clever when I first implemented this, but a
        definite side-effect is the necessity of an overly complicated
        __setattr__() method to make sure that the underlying __data_dict
        class member is properly updated when we change the table values---which
        we do all the time, e.g., when summing binned products.

        Note that, as of version 12, this only works if the attribute name is
        spelled upper-case. This is intented to make more obvious in the source
        code when we are trying to access data through this __getattr__ hook,
        vs. normal member access---I think it was a terrible mistake to make this
        case-insensitive to start with. Note I added some basic logic to give a
        sensible error message when a user try an old-style, lower-case access.
        """
        if not name.isupper():
            logger.error('%s has no attribute %s.', self.__class__.__name__, name)
            if name.upper() in self.__data_dict:
                msg = 'Please use the upper-case form "%s" (new in version 12.0.0).'
                logger.info(msg, name.upper())
            abort()
        return self.__data_dict[name]

    def __setattr__(self, name, value):
        """Proxy method to set an entry in the underlying data as if they were
        class members.

        Warning
        -------
        This is admittedly quite terrible, but I only realize the need for an
        automatic mechanism for updating __data_dict very late in the game, and
        now the data access by name is all over the place.

        A couple of things to remember:

        * within this function we need to access the class member via a call to
          the parent object.__getattribute__() in order to avoid an infinite
          recursion;
        * the __data_dict class member does not always exist in the life of the
          object, and therefore we need to encapsulate the call into a
          try/except block.
        """
        try:
            _data_dict = object.__getattribute__(self, '_xBinnedFileBase__data_dict')
            if name in _data_dict:
                _data_dict[name] = value
                # Update the underlyng HDU list---this is important to have
                # the thing always up to date in case you want to write it
                # back to file at some point. This presents the additional
                # challenge that we have to map back the __data_dict class
                # member to the actual structure in the underlying FITS files.
                # The following hack only works as long as we are dealing with
                # either image extensions or binary tables (at position 1).
                # This might be another reason why we might want to refactor the
                # entire xBinnedFileBase class at some point.
                try:
                    self.hdu_list[name].data = value
                except KeyError:
                    self.hdu_list[1].data[name] = value
        except AttributeError:
            pass
        object.__setattr__(self, name, value)

    def set_data(self, name, value):
        """Add a (key, value) pair to the underlying __data_dict class member.

        Warning
        -------
        This is yet another side effect of the perverse __getattr__ / __setattr__
        mechanism we have put in place for this class---if we calculate a
        quantity dinamically for a binned object and we still want to be able
        to access it with the same rules, we have to manually add it to the
        dictionary.
        """
        assert name.isupper()
        self.__data_dict[name] = value

    def _check_iadd(self, other, same_shape=(), same_values=()):
        """Convenience method to check that two binned object are suitable for
        addition.

        Typically we want to make sure that the two object are of the same class
        and that the underlying data arrays have the same shape. Occasionally,
        we also require that some of the data, e.g., the EBONDS, are identical.

        This is a convenient hook to have for the implementation of the
        __iadd__() dunder of the subclasses.
        """
        assert isinstance(other, self.__class__)
        for key in same_shape:
            s1 = self.__data_dict[key].shape
            s2 = other.__data_dict[key].shape
            if s1 != s2:
                logger.error('Arrays have different shapes (%s vs %s) for key "%s"', s1, s2, key)
                abort('Binned files cannot be combined consistently')
        for key in same_values:
            v1 = self.__data_dict[key]
            v2 = other.__data_dict[key]
            if not numpy.array_equal(v1, v2):
                logger.error('Arrays have different values for key "%s"', key)
                logger.error('First array:\n%s', v1)
                logger.error('Second array:\n%s', v2)
                try:
                    delta = v1 - v2
                    logger.error('Difference:\n%s', delta)
                except ValueError as e:
                    logger.error('Cannot calculate difference: %s', e)
                abort('Binned files cannot be combined consistently')

    def ontime(self, scale=1.e-3):
        """Return the nominal ontime for the binned file.

        Arguments
        ---------
        scale : float
            A scale factor, by default 1.e-3 (i.e., ontime in ks).
        """
        return self.primary_header['ONTIME'] * scale

    def binning_algorithm(self):
        """Return the binning algorithm used to create the file.
        """
        return self.primary_header['BINALG']

    def du_id(self):
        """Return the DU ID for the binned file.

        .. versionadded:: 12.0.0

        Warning
        -------
        This was added in version 12 of ixpeobssim after
        https://bitbucket.org/ixpesw/ixpeobssim/issues/327
        """
        key = 'DETNAM'
        try:
            return int(self.primary_header[key][-1])
        except KeyError:
            logger.error('The binned file primary header has no %s keyword.', key)
            logger.error('(Note we did not propagate %s to binned file prior to version 12)', key)
            return None

    # pylint: disable=unused-argument
    def plot(self, *args, **kwargs):
        """Do nothing method.
        """
        logger.info('%s.plot() not implemented, yet.', self.__class__.__name__)

    def write(self, file_path, overwrite=True):
        """Write the binned data to the output file.
        """
        logger.info('Writing %s object to %s...', self.__class__.__name__, file_path)
        self.hdu_list.writeto(file_path, overwrite=overwrite)
        logger.info('Done.')

    def __str__(self):
        """String formatting.
        """
        return '%s content:\n%s' % (self.__class__.__name__, self.__data_dict)



def peek_binning_algorithm(file_path):
    """Open a binned file and peek at the underlying algorithm
    """
    return xBinnedFileBase(file_path).binning_algorithm()

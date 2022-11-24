# Copyright (C) 2022, the ixpeobssim team.
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

"""facilities for exposure calculation.
"""

from __future__ import print_function, division

from astropy.io import fits
from astropy import wcs
import numpy


from ixpeobssim.binning.base import xEventBinningBase, xBinnedFileBase
from ixpeobssim.binning.fmt import xBinTableHDUEBOUNDS, xBinTableHDUTHETABOUNDS
from ixpeobssim.core.fitsio import xFITSImageBase, xPrimaryHDU
from ixpeobssim.core.hist import xHistogram3d
from ixpeobssim.instrument.gpd import within_fiducial_rectangle
from ixpeobssim.instrument.mma import sky_to_gpd
from ixpeobssim.utils.astro import region_compound, angular_separation
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt

# pylint: disable=invalid-name, too-many-arguments, too-many-locals, no-member


class xEventBinningLTCUBE(xEventBinningBase):

    """Class for LTCUBE binning.
    """

    INTENT = 'Livetime map cube in sky coordinates and off axis angle bins'
    SUPPORTED_KWARGS = xEventBinningBase._image_wcs_kwargs()
    #EXTENSIONS_NAMES = ['LIVETIME', 'THETA_BOUNDS']

    def process_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.process_kwargs(self)
        self.process_image_ref_kwargs()

    def bin_(self):
        """Overloaded method.
        """
        # Retrieve the spacecraft extension from the input event file.
        sc_data = self.event_file.sc_data()
        # Get the target occult column and create a mask to select when the
        # target is visible and spacecraft is not in SAA.
        # NB: We are getting rid of the last item (here and after) because we are
        # treating each MET as start time of the interval.
        mask = numpy.logical_and(sc_data['TARGET_OCCULT'][:-1] == 0,
                                     sc_data['IN_SAA'][:-1] == 0)
        met = sc_data['MET']
        # Here we should use the livetime per interval, not just the delta time.
        livetime = (met[1:] - met[:-1])[mask]
        ra_pnt = sc_data['RA_PNT'][:-1][mask]
        dec_pnt = sc_data['DEC_PNT'][:-1][mask]
        du_id = self.event_file.du_id()
        # Build the wcs for the final maps
        wcs_ = self._build_image_wcs(**self.kwargs)
        nx, ny = wcs_.array_shape
        # Create the binning for the 3d histogram (off-axis angle, x, y).
        xbinning = numpy.linspace(0, nx, nx + 1).astype(int)
        ybinning = numpy.linspace(0, ny, ny + 1).astype(int)
        thetabins = self.kwargs.get('thetabins')
        zbinning = numpy.linspace(0, 8.5, thetabins + 1)
        histo3d = xHistogram3d(zbinning, xbinning, ybinning)
        # Build a meshgrid of pixel indexes from the bin edges (indexing=xy or ij?)
        x, y = numpy.meshgrid(xbinning[:-1], ybinning[:-1], indexing='xy')
        x = x.flatten()
        y = y.flatten()
        # Transform the pixel centers to sky coordinates using wcs.
        # TODO: not sure why we need to swap x and y here to get a DU rotation
        # coherent with the CMAP. Are we doing things right?
        ra, dec = wcs_.pixel_to_world_values(y, x)
        # Now we loop over the item of SC_DATA and fill the histogram.
        logger.info('Creating LTCUBE...')
        for lvt, _ra_pnt, _dec_pnt in zip(livetime, ra_pnt, dec_pnt):
            # Calculate the angular separation for all the pixels of the map
            theta = angular_separation(ra, dec, _ra_pnt, _dec_pnt) * 60. #arcmin
            # Find the indexes for the off-axis part of the histogram. We remove
            # 1 so that the first bin is 0.
            z = numpy.digitize(theta, zbinning) - 1
            # We should pass the roll angle here
            detx, dety = sky_to_gpd(ra, dec, None, _ra_pnt, _dec_pnt, du_id, 0.)
            # Create a mask cutting the GPD fiducial area and the angles greater
            # than the uppermost edge of the theta binning (for which
            # z==len(zbinning)).
            mask = (z != thetabins) * within_fiducial_rectangle(detx, dety)
            histo3d.content[z[mask], x[mask], y[mask]] += lvt
        # Write all the info to output fits file.
        data = histo3d.content
        header = wcs_.to_header()
        hdu_list = fits.HDUList()
        hdu = self.build_primary_hdu(data, header=header)
        hdu.header['EXTNAME'] = 'LIVETIME'
        hdu_list.append(hdu)
        logger.info('Creating THETA BOUNDS...')
        zbounds = xBinTableHDUTHETABOUNDS((zbinning[:-1], zbinning[1:]))
        zbounds.setup_header(self.event_file.primary_keywords())
        hdu_list.append(zbounds)
        self.write_output_file(hdu_list)



class xBinnedLivetimeCube(xBinnedFileBase):

    """ Read-mode interface to a LTCUBE FITS file.

    TODO: Temporarily most of this class is just a copy of xBinnedMDPMapCube.
          We should refactor the common methods in a separate class.
    """

    def _read_data(self):
        """Overloaded method.
        """
        self._img_header = self.hdu_list['LIVETIME'].header
        self.wcs = wcs.WCS(self._img_header)
        self._read_image_data('LIVETIME')
        self._read_binary_table_data('THETA_BOUNDS')

    def __iadd__(self, other):
        """Overloaded method for binned data addition.
        """
        same_shape = ['LIVETIME']
        same_values = ('THETA_LO','THETA_HI')
        self._check_iadd(other, same_shape, same_values)
        self.LIVETIME += other.LIVETIME
        return self

    def map_shape(self):
        """Return the shape of the underlying sky-maps.

        Mind the underlying arrays are all 3-dimensional cubes with the theta
        binning as the first dimension---so it's the last two that we care about.
        (And, for completeness: since the arrays are all the same, we use I
        for convenience.)
        """
        print('shape = ',self.LIVETIME.shape[1:])
        return self.LIVETIME.shape[1:]

    def pixel_size(self):
        """Return the pixel size of the underlying spatial map.
        """
        cdelt1, cdelt2, _ = self.wcs.wcs.get_cdelt()
        assert abs(cdelt1) == abs(cdelt2)
        return abs(cdelt1)

    def ds9_region_mask(self, region_list, region_slice=None):
        """Return the (spatial) array map corresponding to a given ds9 region list.

        If there is more than one region in the region list, by default the
        output mask corresponds to the logical or of all the regions.
        The region_slice optional argument allows to restrict the mask to a
        subset of the regions.

        Arguments
        ---------
        region_list :
            The region list defining the mask

        region_slice : int or slice (optional)
            An optional slice designator to select a subset of the region list.
        """
        if region_slice is not None:
            if isinstance(region_slice, int):
                region_slice = slice(region_slice, region_slice + 1)
            region_list = region_list[region_slice]

        region = region_compound(region_list)[0].to_pixel(self.wcs)

        return region.to_mask().to_image(self.map_shape()).astype(bool)

    @staticmethod
    def _sum_image_pixels(data, spatial_mask, theta_layer=0):
        """Sum one of the underlying image extentions over a given spatial
        mask and theta slice.
        """
        return data[theta_layer][spatial_mask].sum()

    def sum_pixels(self, spatial_mask, theta_layer=0):
        """Sum the relevant quantities over a given spatial mask and theta slice.
        """
        args = spatial_mask, theta_layer
        LIVETIME = self._sum_image_pixels(self.LIVETIME, *args)
        return {'LIVETIME': LIVETIME}

    def theta_centers(self):
        """ Return the centers of the binning over theta.
        """
        return 0.5 * (self.THETA_HI + self.THETA_LO)

    def theta_binning(self):
        """Return the underlying theta binning in the form of a numpy array.

        Note the array has one more element than the underlying THETA_LO and
        THETA_HI arrays, i.e., if the cube is binned in two theta layers,
        say 0--2 arcmin and 2--4 arcmin, the theta binning is [0. 2. 4.].
        """
        return numpy.append(self.THETA_LO, self.THETA_HI[-1])

    def num_theta_layers(self):
        """Return the number of theta layers in the cube.
        """
        return len(self.THETA_LO)

    def theta_range(self, theta_layer):
        """Return the theta value for a given theta layer.
        """
        return self.THETA_LO[theta_layer],self.THETA_HI[theta_layer]

    def theta_label(self, theta_layer):
        """Return a text label corresponding to the theta layer (e.g, to
        identify the theta range on a plot).
        """
        return '%.2f--%.2f arcmin' % self.theta_range(theta_layer)

    def _plot_layer(self, data, theta_layer, **kwargs):
        """Delegated function to make a plot of one of the underlying image
        extensions in *one* theta layer of the cube.

        This essentially calculates the proper slice of the input data and all
        the necessary arguments to be passed to xFITSImageBase for doing the
        actual plot.

        Note this function is *not* creating any matplotlib canvas---it is left
        to the caller to do that.

        Arguments
        ---------
        data : array
            The underlying data array. This can be the array corresponding to
            one of the underlying image extensions of the cube, or any
            array of the same shape.

        theta_range : int
            The identifier of the theta layer in the cube.

        kwargs : dict
            Additional keyword arguments being passed to the underlying
            xFITSImageBase.make_plot() call.
        """
        assert theta_layer in range(0, self.num_theta_layers())
        data = data[theta_layer, :, :]
        slices = ('x', 'y', theta_layer)
        xFITSImageBase.make_plot(data, self.wcs, slices=slices, **kwargs)
        xFITSImageBase.add_label(self.theta_label(theta_layer), y=0.92)

    def _plot_base(self, data, theta_layers, mask, figname, prefix=None,
                   zlabel=None, post_plot_hook=None, **kwargs):
        """Delegated base function to plot multiple theta layers of the same
        image extension.
        """
        # Do some magic with the theta_layers argument.
        if theta_layers is None:
            theta_layers = list(range(0, self.num_theta_layers()))
        elif isinstance(theta_layers, int):
            theta_layers = [theta_layers]
        # If we are passing a mask, we create a copy of the data array (not to)
        # mess up with the underlying data, and explicitely set to zero all the
        # values that are not in the mask.
        if mask is not None:
            data = data.copy()
            data[numpy.logical_not(mask)] = 0.
        # Loop over the theta layers.
        if zlabel is None:
            zlabel = figname
        kwargs.setdefault('zlabel', zlabel)
        for theta_layer in theta_layers:
            figure_name = '%s %s' % (figname, self.theta_label(theta_layer))
            if prefix is not None:
                figure_name = '%s %s' % (prefix, figure_name)
            plt.figure(figure_name)
            self._plot_layer(data, theta_layer, **kwargs)
            if post_plot_hook is not None:
                post_plot_hook(theta_layer, mask)

    def plot_livetime_map(self, theta_layers=None, prefix=None, **kwargs):
        """Plot the MDP map.
        """
        self._plot_base(self.LIVETIME, theta_layers, None, 'LIVETIME', prefix, 'LIVETIME', **kwargs)

    def plot(self, prefix=None):
        """Plot the data.
        """
        self.plot_livetime_map(prefix=prefix)



class xExposureCube:

    """Structure for writing/reading exposure cubes.
    """

    def __init__(self, exposure, ebounds=None, units=None):
        self.exposure = exposure
        self.units = units
        if ebounds is not None:
            self.energy_lo, self.energy_hi = ebounds
            self.energy = 0.5 * (self.energy_lo + self.energy_hi)

    @classmethod
    def empty(cls, shape):
        """Create an empty exposure cube.
        """
        exposure = numpy.zeros(shape)
        return cls(exposure)

    @classmethod
    def from_file(cls, file_path):
        """Read an exposure cube from file.
        """
        raise NotImplementedError

    def write(self, file_path, header, overwrite):
        """Write the exposure cube to file.

        .. warning::

           Need to add all the keywords from the LTCUBE file
        """
        hdu_list = fits.HDUList()
        hdu = xPrimaryHDU(self.exposure, header, 'xpexposure.py')
        hdu.set_keyword('EXTNAME', 'EXPOSURE')
        hdu.add_keyword('BUNIT', self.units, 'Units of image array values')
        hdu_list.append(hdu)
        logger.info('Creating ENERGY BOUNDS...')
        ebounds = xBinTableHDUEBOUNDS((self.energy_lo, self.energy_hi))
        hdu_list.append(ebounds)
        hdu_list.writeto(file_path, overwrite=overwrite)



def create_psf_kernel(psf, pixsize, npix=21):
    """Generate a binned kernel for the PSF convolution.

    .. warning::

       There is potential overlap with the calculate_psf_kernel() function in
       the evt.deconvolution module, but the latter is doing something
       slightly different---a Monte Carlo integration over the pixels, rather than
       a direct evaluation at the pixel center. We should probably think about
       whether which one we want.
    """
    psf_kernel = numpy.zeros((1, npix, npix))
    for i in range(npix):
        for j in range(npix):
            x = (i - npix/2. + 0.5) * pixsize * 3600
            y = (j - npix/2. + 0.5) * pixsize * 3600
            d = numpy.sqrt(x**2 + y**2)
            psf_kernel[0,i,j] = psf(d)
    psf_kernel /= psf_kernel.sum()
    return psf_kernel

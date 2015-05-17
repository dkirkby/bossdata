# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

""" Access spectroscopic data for a single BOSS target.
"""

from __future__ import division,print_function

import numpy as np
import numpy.ma

import fitsio

from bossdata.bits import SPPIXMASK

class SpecFile(object):
    """ A BOSS spec file containing summary data for a single target.

    The spec file data model is described at
    http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html.
    This class supports the full version described in the data model as well as a "lite"
    version that does not contain the per-exposure HDUs with indices >= 4. Use the `lite`
    attribute to detect which version an object represents.

    This class is only intended for reading the BOSS spec file format, so generic
    operations on spectroscopic data (redshifting, etc) are intentionally not included here.

    Args:
        path(str): Local path of the spec FITS file to use.  This should normally be obtained
            via :meth:`bossdata.path.Finder.get_spec_path` and can be automatically mirrored via
            :meth:`bossdata.remote.Manager.get` or using the :ref:`bossfetch` script. The file
            is opened in read-only mode so you do not need write privileges.
    """
    def __init__(self,path):
        self.hdulist = fitsio.FITS(path,mode = fitsio.READONLY)
        self.lite = (len(self.hdulist) == 4)

    def get_spectrum_hdu(self,exposure_index=None):
        """Get the HDU containing a specified spectrum.

        Args:
            exposure_index(int): Individual exposure to use. Uses the co-added spectrum when
                this is None.

        Returns:
            The HDU object corresponding to the requested spectrum.

        Raises:
            ValueError: Invalid exposure_index.
        """
        if exposure_index is None:
            return self.hdulist[1]
        if exposure_index < 0:
            raise ValueError('exposure index must be >= 0.')
        try:
            return self.hdulist[4+exposure_index]
        except IndexError:
            raise ValueError('Invalid exposure index {0:d}. Valid range is 0-{1:d}.'.format(
                exposure_index,len(self.hdulist)-4))

    def get_valid_data(self,exposure_index=None,pixel_quality_mask=None):
        """Get the valid for a specified exposure or the combined coadd.

        You will probably find yourself using this idiom often::

            data = spec.get_valid_data(...)
            wlen,flux,dflux = data['wavelength'][:],data['flux'][:],data['dflux'][:]

        Args:
            exposure_index(int): Individual exposure to use. Uses the co-added spectrum when
                this is None.
            pixel_quality_mask(int): An integer value interpreted as a bit pattern using the
                bits defined in :attr:`bossdata.bits.SPPIXMASK` (see also
                http://www.sdss3.org/dr10/algorithms/bitmask_sppixmask.php). Any bits set in
                this mask are considered harmless and the corresponding spectrum pixels are
                assumed to contain valid data. When accessing the coadded spectrum, this mask
                is applied to the AND of the masks for each individual exposure. No mask is
                applied if this value is None.

        Returns:
            numpy.ma.MaskedArray: Masked array of per-pixel records. Pixels with no valid data
                are included but masked. The record for each pixel has three named fields:
                wavelength in Angstroms, flux and dflux in 1e-17 ergs/s/cm2/Angstrom. Wavelength
                values are strictly increasing and dflux is calculated as ivar**-0.5 for pixels
                with valid data.
        """
        hdu = self.get_spectrum_hdu(exposure_index)

        if exposure_index is None:
            pixel_bits = hdu['and_mask'][:]
        else:
            pixel_bits = hdu['mask'][:]
        if pixel_quality_mask is not None:
            clear_allowed = np.bitwise_not(np.uint32(pixel_quality_mask))
            pixel_bits = np.bitwise_and(pixel_bits,clear_allowed)
        num_pixels = len(pixel_bits)

        # Identify the pixels with valid data.
        ivar = hdu['ivar'][:]
        bad_pixels = (pixel_bits != 0) | (ivar <= 0.0)
        good_pixels = ~bad_pixels

        # Create and fill the unmasked structured array of data.
        data = np.empty(num_pixels,dtype = [
            ('wavelength',np.float32),('flux',np.float32),('dflux',np.float32)])
        data['wavelength'][:] = np.power(10.0,hdu['loglam'][:])
        data['flux'][:] = hdu['flux'][:]
        data['dflux'][good_pixels] = 1.0/np.sqrt(ivar[good_pixels])
        data['dflux'][bad_pixels] = 0.0

        return np.ma.MaskedArray(data,mask=bad_pixels)

# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

""" Access spectroscopic data for a single BOSS target.
"""

from __future__ import division, print_function

import re

import numpy as np
import numpy.ma

import fitsio

import astropy.table


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
            via :meth:`bossdata.path.Finder.get_spec_path` and can be automatically mirrored
            via :meth:`bossdata.remote.Manager.get` or using the :ref:`bossfetch` script. The
            file is opened in read-only mode so you do not need write privileges.
    """
    def __init__(self, path):
        self.hdulist = fitsio.FITS(path, mode=fitsio.READONLY)
        self.lite = (len(self.hdulist) == 4)
        self.header = self.hdulist[0].read_header()
        # Look up the available exposures.
        self.num_exposures = self.header['NEXP']
        exposures = {}
        expid_pattern = re.compile('([br][12])-([0-9]{8})-([0-9]{8})-([0-9]{8})')
        for i in range(self.num_exposures):
            spec_id, exp_num, flat_num, arc_num = expid_pattern.match(
                self.header['EXPID{0:02d}'.format(i + 1)]).groups()
            if exp_num in exposures:
                info = exposures[exp_num]
                # Check that the arc and flat exposure numbers agree.
                if info['arc'] != arc_num:
                    raise RuntimeError(
                        'Found different red/blue arcs for expid {}.'.format(exp_num))
                if info['flat'] != flat_num:
                    raise RuntimeError(
                        'Found different red/blue flats for expid {}.'.format(exp_num))
            else:
                # Initialize a new record with zeros to indicate a missing camera.
                info = dict(arc=arc_num, flat=flat_num, bhdu=0, rhdu=0)
            # Record the HDU number for this camera.
            info[spec_id[0]+'hdu'] = 4 + i
            exposures[exp_num] = info
        # Build a table of exposure info sorted by exposure number.
        self.exposure_table = astropy.table.Table(
            names=('exp', 'arc', 'flat', 'bhdu', 'rhdu'),
            dtype=('i4', 'i4', 'i4', 'i4', 'i4'))
        for exp_num in sorted(exposures.keys()):
            info = exposures[exp_num]
            self.exposure_table.add_row(
                (exp_num, info['arc'], info['flat'], info['bhdu'], info['rhdu']))

    def get_exposure_hdu(self, exposure_index=None, camera=None):
        """Lookup the HDU for one exposure.

        Args:
            exposure_index(int): Individual exposure to use, specified as a sequence number
                starting from zero, for the first exposure, and increasing up to
                `self.num_exposures-1`.
            camera(str): Which camera to use. Must be either 'blue' or 'red'.

        Returns:
            hdu: The HDU containing data for the requested exposure.

        Raises:
            ValueError: Invalid exposure_index or camera.
            RuntimeError: This camera is missing for this exposure.
        """
        if exposure_index < 0 or exposure_index >= self.num_exposures:
            raise ValueError('exposure index must be in the range 0-{0}'.format(
                self.num_exposures - 1))
        if camera not in ('blue', 'red'):
            raise ValueError('camera must be either "blue" or "red".')

        info = self.exposure_table[exposure_index]
        if camera == 'blue':
            hdu_index = info['bhdu']
        else:
            hdu_index = info['rhdu']
        if hdu_index == 0:
            raise RuntimeError('Missing {0} camera for exposure {1}.'.format(
                camera, info['exp']))

        return self.hdulist[hdu_index]

    def get_pixel_mask(self, exposure_index=None, camera=None):
        """Get the pixel mask for a specified exposure or the combined coadd.

        Returns the `and_mask` for coadded spectra. The entire mask is returned, including
        any pixels with zero inverse variance.

        Args:
            exposure_index(int): Individual exposure to use, specified as a sequence number
                starting from zero, for the first exposure, and increasing up to
                `self.num_exposures-1`. Uses the co-added spectrum when the value is None.
            camera(str): Which camera to use. Must be either 'b' (blue) or 'r' (red) unless
                exposure_index is None, in which case this argument is ignored.

        Returns:
            numpy.ndarray: Array of integers, one per pixel, encoding the mask bits defined
                in :attr:`bossdata.bits.SPPIXMASK` (see also
                http://www.sdss3.org/dr10/algorithms/bitmask_sppixmask.php).
        """
        if exposure_index is None:
            hdu = self.hdulist[1]
            return hdu['and_mask'][:]
        else:
            hdu = self.get_exposure_hdu(exposure_index, camera)
            return hdu['mask'][:]

    def get_valid_data(self, exposure_index=None, camera=None, pixel_quality_mask=None,
                       include_wdisp=False, include_sky=False):
        """Get the valid data for a specified exposure or the combined coadd.

        You will probably find yourself using this idiom often::

            data = spec.get_valid_data(...)
            wlen,flux,dflux = data['wavelength'][:],data['flux'][:],data['dflux'][:]

        Args:
            exposure_index(int): Individual exposure to use, specified as a sequence number
                starting from zero, for the first exposure, and increasing up to
                `self.num_exposures-1`. Uses the co-added spectrum when the value is None.
            camera(str): Which camera to use. Must be either 'b' (blue) or 'r' (red) unless
                exposure_index is None, in which case this argument is ignored.
            pixel_quality_mask(int): An integer value interpreted as a bit pattern using the
                bits defined in :attr:`bossdata.bits.SPPIXMASK` (see also
                http://www.sdss3.org/dr10/algorithms/bitmask_sppixmask.php). Any bits set in
                this mask are considered harmless and the corresponding spectrum pixels are
                assumed to contain valid data. When accessing the coadded spectrum, this mask
                is applied to the AND of the masks for each individual exposure. No mask is
                applied if this value is None.
            include_wdisp: Include a wavelength dispersion column in the returned data.
            include_sky: Include a sky flux column in the returned data.

        Returns:
            numpy.ma.MaskedArray: Masked array of per-pixel records. Pixels with no valid data
                are included but masked. The record for each pixel has at least the following
                named fields: wavelength in Angstroms, flux and dflux in 1e-17
                ergs/s/cm2/Angstrom. Wavelength values are strictly increasing and dflux is
                calculated as ivar**-0.5 for pixels with valid data. Optional fields are
                wdisp in Angstroms and sky in 1e-17 ergs/s/cm2/Angstrom.
        """
        # Look up the HDU for this spectrum and its pixel quality bitmap.
        if exposure_index is None:
            hdu = self.hdulist[1]
            pixel_bits = hdu['and_mask'][:]
        else:
            hdu = self.get_exposure_hdu(exposure_index, camera)
            pixel_bits = hdu['mask'][:]
        num_pixels = len(pixel_bits)

        # Apply the pixel quality mask, if any.
        if pixel_quality_mask is not None:
            clear_allowed = np.bitwise_not(np.uint32(pixel_quality_mask))
            pixel_bits = np.bitwise_and(pixel_bits, clear_allowed)

        # Identify the pixels with valid data.
        ivar = hdu['ivar'][:]
        bad_pixels = (pixel_bits != 0) | (ivar <= 0.0)
        good_pixels = ~bad_pixels

        # Create and fill the unmasked structured array of data.
        dtype = [('wavelength', np.float32), ('flux', np.float32), ('dflux', np.float32)]
        if include_wdisp:
            dtype.append(('wdisp', np.float32))
        if include_sky:
            dtype.append(('sky', np.float32))
        data = np.empty(num_pixels, dtype=dtype)
        data['wavelength'][:] = np.power(10.0, hdu['loglam'][:])
        data['flux'][:] = hdu['flux'][:]
        data['dflux'][good_pixels] = 1.0 / np.sqrt(ivar[good_pixels])
        data['dflux'][bad_pixels] = 0.0
        if include_wdisp:
            data['wdisp'][:] = np.power(10., hdu['wdisp'][:])
        if include_sky:
            data['sky'][:] = hdu['sky'][:]

        return numpy.ma.MaskedArray(data, mask=bad_pixels)

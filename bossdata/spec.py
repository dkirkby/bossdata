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


def get_fiducial_pixel_index(wavelength):
        """
        Convert a wavelength to a fiducial pixel index.

        Args:
            wavelength(float): Input wavelength in Angstroms.

        Returns:
            numpy.ndarray: Array of floating-point indices relative to the fiducial
                wavelength grid.
        """
        return (np.log10(wavelength) - _fiducial_log10lam0)/_fiducial_coef

_fiducial_coef = 1e-4
_fiducial_log10lam0 = np.log10(3500.26)

fiducial_pixel_index_range = (0, 4800)
"""
Range of fiducial pixel indices covering all spectra.

Pixel indices are calculated by :func:`get_fiducial_pixel_index`.
"""

fiducial_loglam = (_fiducial_log10lam0 +
                   _fiducial_coef * np.arange(*fiducial_pixel_index_range))
"""
Array of fiducial log10(wavelength in Angstroms) covering all spectra.
"""


class Exposures(object):
    """Table of exposure info extracted from FITS header keywords.

    Uses the NEXP and EXPIDnn keywords that are present in the header of HDU0
    in spPlate and spec FITS files.

    Args:
        header(dict): dictionary of FITS header keyword, value pairs.

    Returns:
        `astropy.table.Table`: table with columns 'exp', 'arc', 'flat', 'b1', 'b2',
        'r1', 'r2', that give the science exposure numbers 'exp', together with the
        corresponding 'arc' and 'flat' exposure numbers used to calibrate them.
        The 'b1', 'b2', 'r1', 'r2' columns give the value n=1-12 of the 'EXP<n>' header
        keyword described each of the four spectrographs, or n=0 if this spectrograph
        is not present for this exposure.
    """
    def __init__(self, header):
        num_exposures = header['NEXP']
        expid_pattern = re.compile('([br][12])-([0-9]{8})-([0-9]{8})-([0-9]{8})')
        exposure_set = set()
        self.table = astropy.table.Table(
            names=('offset', 'camera', 'science', 'flat', 'arc'),
            dtype=('i4', 'S2', 'i4', 'i4', 'i4'))
        for i in range(num_exposures):
            camera, science_num, flat_num, arc_num = expid_pattern.match(
                header['EXPID{0:02d}'.format(i + 1)]).groups()
            self.table.add_row((i, camera, int(science_num), int(flat_num), int(arc_num)))
            exposure_set.add(int(science_num))
        self.sequence = sorted(exposure_set)
        self.num_exposures = len(exposure_set)

    def get_info(self, exposure_index, camera):
        """Get info about the specified exposure.

        Args:
            exposure_index(int): The sequence number for the requested exposure,
                starting from zero.
            camera(str): One of b1,b2,r1,r2.

        Returns:
            A structured array with information about the requested exposure.

        Raises:
            ValueError: Invalid exposure_index or camera.
            RuntimeError: Exposure not present.
        """
        if camera not in ('b1', 'b2', 'r1', 'r2'):
            raise ValueError(
                'Invalid camera "{}", expected b1, b2, r1, or r2.'.format(camera))
        if exposure_index < 0 or exposure_index >= self.num_exposures:
            raise ValueError('Invalid exposure_index {}, expected 0-{}.'.format(
                exposure_index, self.num_exposures))
        science_num = self.sequence[exposure_index]
        row = (self.table['science'] == science_num) & (self.table['camera'] == camera)
        if not np.any(row):
            raise RuntimeError('No exposure[{}] = {:08d} found for {}.'.format(
                exposure_index, science_num, camera))
        if np.count_nonzero(row) > 1:
            raise RuntimeError(
                'Found multiple {} exposures {:08d}.'.format(camera, exposure_index))
        return self.table[row][0]


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
        self.exposures = Exposures(self.header)

    def get_exposure_hdu(self, exposure_index, camera):
        """Lookup the HDU for one exposure.

        This method will not work on "lite" files, which do not include individual
        exposures.

        Args:
            exposure_index(int): Individual exposure to use, specified as a sequence number
                starting from zero, for the first exposure, and increasing up to
                `self.num_exposures-1`.
            camera(str): Which camera to use. Must be one of b1,b2,r1,r2.

        Returns:
            hdu: The HDU containing data for the requested exposure.

        Raises:
            RuntimeError: individual exposures not available in lite file.
        """
        if self.lite:
            raise RuntimeError('individual exposures not available in lite file.')

        info = self.exposures.get_info(exposure_index, camera)
        return self.hdulist[4 + info['offset']]

    def get_pixel_mask(self, exposure_index=None, camera=None):
        """Get the pixel mask for a specified exposure or the combined coadd.

        Returns the `and_mask` for coadded spectra. The entire mask is returned, including
        any pixels with zero inverse variance.

        Args:
            exposure_index(int): Individual exposure to use, specified as a sequence number
                starting from zero, for the first exposure, and increasing up to
                `self.num_exposures-1`. Uses the co-added spectrum when the value is None.
            camera(str): Which camera to use. Must be either 'b1', 'b2' (blue) or 'r1', 'r2'
                (red) unless exposure_index is None, in which case this argument is ignored.

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
                       include_wdisp=False, include_sky=False, use_ivar=False,
                       use_loglam=False, fiducial_grid=False):
        """Get the valid data for a specified exposure or the combined coadd.

        You will probably find yourself using this idiom often::

            data = spec.get_valid_data(...)
            wlen,flux,dflux = data['wavelength'][:],data['flux'][:],data['dflux'][:]

        Args:
            exposure_index(int): Individual exposure to use, specified as a sequence number
                starting from zero, for the first exposure, and increasing up to
                `self.num_exposures-1`. Uses the co-added spectrum when the value is None.
            camera(str): Which camera to use. Must be either 'b1', 'b2' (blue) or 'r1', 'r2'
                (red) unless exposure_index is None, in which case this argument is ignored.
            pixel_quality_mask(int): An integer value interpreted as a bit pattern using the
                bits defined in :attr:`bossdata.bits.SPPIXMASK` (see also
                http://www.sdss3.org/dr10/algorithms/bitmask_sppixmask.php). Any bits set in
                this mask are considered harmless and the corresponding spectrum pixels are
                assumed to contain valid data. When accessing the coadded spectrum, this mask
                is applied to the AND of the masks for each individual exposure. No mask is
                applied if this value is None.
            include_wdisp: Include a wavelength dispersion column in the returned data.
            include_sky: Include a sky flux column in the returned data.
            use_ivar: Replace ``dflux`` with ``ivar`` (inverse variance) in the returned
                data.
            use_loglam: Replace ``wavelength`` with ``loglam`` (``log10(wavelength)``) in
                the returned data.
            fiducial_grid: Return co-added data using the fiducial wavelength grid defined
                by :attr:`fiducial_pixel_index_range`.  If False, the returned array uses
                the native grid of the SpecFile, which generally trims pixels on both ends
                that have zero inverse variance.  Set this value True to ensure that all
                co-added spectra use aligned wavelength grids when this matters.

        Returns:
            numpy.ma.MaskedArray: Masked array of per-pixel records. Pixels with no valid data
                are included but masked. The record for each pixel has at least the following
                named fields: wavelength in Angstroms (or loglam), flux and dflux in 1e-17
                ergs/s/cm2/Angstrom (or flux and ivar). Wavelength values are strictly
                increasing and dflux is calculated as ivar**-0.5 for pixels with valid data.
                Optional fields are wdisp in constant-log10-lambda pixels and sky in 1e-17
                ergs/s/cm2/Angstrom.
        """
        # Look up the HDU for this spectrum and its pixel quality bitmap.
        if exposure_index is None:
            hdu = self.hdulist[1]
            pixel_bits = hdu['and_mask'][:]
        else:
            hdu = self.get_exposure_hdu(exposure_index, camera)
            pixel_bits = hdu['mask'][:]

        if fiducial_grid:
            if exposure_index is not None:
                raise ValueError('Fiducial grid not supported for individual exposures.')
            loglam = fiducial_loglam
            first_index = float(get_fiducial_pixel_index(10.0**hdu['loglam'][0]))
            if abs(first_index - round(first_index)) > 0.01:
                raise RuntimeError('Wavelength grid not aligned with fiducial grid.')
            trimmed = slice(first_index, first_index + len(pixel_bits))
        else:
            loglam = hdu['loglam'][:]
            trimmed = slice(None)
        num_pixels = len(loglam)

        # Apply the pixel quality mask, if any.
        if pixel_quality_mask is not None:
            clear_allowed = np.bitwise_not(np.uint32(pixel_quality_mask))
            pixel_bits = np.bitwise_and(pixel_bits, clear_allowed)

        # Identify the pixels with valid data.
        ivar = hdu['ivar'][:]
        bad_pixels = (pixel_bits != 0) | (ivar <= 0.0)
        good_pixels = ~bad_pixels

        # Create and fill the unmasked structured array of data.
        dtype = [('loglam' if use_loglam else 'wavelength', np.float32),
                 ('flux', np.float32), ('ivar' if use_ivar else 'dflux', np.float32)]
        if include_wdisp:
            dtype.append(('wdisp', np.float32))
        if include_sky:
            dtype.append(('sky', np.float32))
        data = np.zeros(num_pixels, dtype=dtype)
        if use_loglam:
            data['loglam'][:] = loglam
        else:
            data['wavelength'][:] = np.power(10.0, loglam)
        data['flux'][trimmed][:] = hdu['flux'][:]
        if use_ivar:
            data['ivar'][trimmed][good_pixels] = ivar[good_pixels]
        else:
            data['dflux'][trimmed][good_pixels] = 1.0 / np.sqrt(ivar[good_pixels])
        if include_wdisp:
            data['wdisp'][trimmed] = hdu['wdisp'][:]
        if include_sky:
            data['sky'][trimmed] = hdu['sky'][:]

        if fiducial_grid:
            mask = np.ones(num_pixels, dtype=bool)
            mask[trimmed][:] = bad_pixels
        else:
            mask = bad_pixels

        return numpy.ma.MaskedArray(data, mask=mask)

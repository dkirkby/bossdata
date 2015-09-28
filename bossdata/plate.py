# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Access BOSS plate data products.
"""

from __future__ import division, print_function

import os.path

import numpy as np
import numpy.ma
import numpy.polynomial.legendre

import astropy.table

import fitsio

from bossdata.spec import Exposures, fiducial_loglam, get_fiducial_pixel_index


def get_num_fibers(plate):
    """Return the number of fiber holes for a given plate number.

    Plate numbers 3510 or larger are (e)BOSS plates with 1000 fibers. Smaller plate
    numbers are assumed to be SDSS-I/II with 640 fibers.

    Args:
        plate(int): Plate number.

    Returns:
        int: The value 640 or 1000.
    """
    return 640 if plate < 3510 else 1000


class Plan(object):
    """The plan file for configuring the BOSS pipeline to combine exposures of a single plate.

    Combined :datamodel:`plan files <PLATE4/spPlan>` are small text files that list the
    per-spectrograph (b1,b2,r1,r2) exposures used as input to a single coadd.  Use the
    `exposure_table` attribute to access this information.  Note that
    :class:`bossdata.spec.SpecFile` has a similar `exposures` attribute which only includes
    exposures actually used in the final co-add, so is generally a subset of the planned
    exposures.

    Args:
        path(str): The local path to a plan file.
    """
    def __init__(self, path):
        self.plate = None
        self.exposures = {}
        with open(path, 'r') as f:
            for line in f:
                if not line.startswith('SPEXP'):
                    continue
                tokens = line.split()
                # token[1] is the plate number and must be identical for all lines.
                if not self.plate:
                    self.plate = int(tokens[1])
                elif self.plate != int(tokens[1]):
                    raise RuntimeError('Internal error: unexpected plate {0} in {1}.'.format(
                        tokens[1], path))
                # Validate the exposure names, which should all have the form
                # [prefix]-[cc]-[eeeeeeee].fits
                exposure_id = set()
                specs = set()
                for name in tokens[7:11]:
                    if name == 'UNKNOWN':
                        continue
                    name, ext = os.path.splitext(name)
                    if ext != '.fits':
                        raise RuntimeError('Unexpected extension {} in {}.'.format(ext, path))
                    fields = name.split('-')
                    if len(fields) != 3:
                        raise RuntimeError(
                            'Unexpected exposure name {} in {}.'.format(name, path))
                    if fields[0] not in ('sdR', 'spFrame'):
                        raise RuntimeError(
                            'Unexpected prefix {} in {}.'.format(fields[0], path))
                    if fields[1] not in ('r1', 'b1', 'r2', 'b2'):
                        raise RuntimeError(
                            'Unexpected camera {} in {}.'.format(fields[1], path))
                    exposure_id.add(int(fields[2]))
                    specs.add(fields[1])
                if len(exposure_id) != 1:
                    raise RuntimeError('Multiple exposure IDs: {}.'.format(exposure_id))
                # Build an exposure record to save.
                exposure = dict(
                    MJD=tokens[2], EXPTIME=float(tokens[5]),
                    EXPID=exposure_id.pop(), SPECS=specs)
                # Record this exposure under the appropriate category.
                flavor = tokens[4]
                if flavor in self.exposures:
                    self.exposures[flavor].append(exposure)
                else:
                    self.exposures[flavor] = [exposure]
        if 'science' in self.exposures:
            self.num_science_exposures = len(self.exposures['science'])
        else:
            self.num_science_exposures = 0
        # Remember the number of fibers on this plate.
        self.num_fibers = get_num_fibers(self.plate)
        # Build a table of science exposure info.
        self.exposure_table = astropy.table.Table(
            names=('exp', 'mjd', 'exptime', 'cameras'), dtype=('i4', 'i4', 'f4', 'S12'))
        for exposure in self.exposures['science']:
            self.exposure_table.add_row((
                exposure['EXPID'], exposure['MJD'], exposure['EXPTIME'],
                ','.join(sorted(exposure['SPECS']))))

    def get_spectrograph_index(self, fiber):
        """Get the spectrograph index 1,2 for the specified fiber.

        Args:
            fiber(int): Fiber number to identify which spectrograph to use, which must
                be in the range 1-1000 (or 1-640 for plate < 3510).

        Returns:
            int: Value of 1 if fiber is read out by the first spectrograph 1-500 (1-320),
                or else 2 for the second spectrograph.

        Raises:
            ValueError: fiber is outside the allowed range 1-1000 (1-640) for this plate.
        """
        if fiber < 1 or fiber > self.num_fibers:
            raise ValueError('Invalid fiber {} should be in the range 1-{} for plate {}.'
                             .format(fiber, self.num_fibers, self.plate))
        return 1 if fiber <= self.num_fibers // 2 else 2

    def get_exposure_name(self, sequence_number, band, fiber, ftype='spCFrame'):
        """Get the file name of a single science exposure data product.

        Use the exposure name to locate FITS data files associated with
        individual exposures.  The supported file types are:
        :datamodel:`spCFrame <PLATE4/spCFrame>`,
        :datamodel:`spFrame <PLATE4/spFrame>`,
        :datamodel:`spFluxcalib <PLATE4/spFluxcalib>` and
        :datamodel:`spFluxcorr <PLATE4/spFluxcorr>`.
        Note that this method returns None when the requested exposure is not present
        in the plan, so the return value should always be checked.

        Args:
            sequence_number(int): Science exposure sequence number, counting from zero.
                Must be less than our num_science_exposures attribute.
            fiber(int): Fiber number to identify which spectrograph to use, which must
                be in the range 1-1000 (or 1-640 for plate < 3510).
            band(str): Must be 'blue' or 'red'.
            ftype(str): Type of exposure file whose name to return.  Must be one of
                spCFrame, spFrame, spFluxcalib, spFluxcorr.  An spCFrame is assumed
                to be uncompressed, and all other files are assumed to be compressed.

        Returns:
            str: Exposure name of the form [ftype]-[cc]-[eeeeeeee].[ext] where [cc]
                identifies the spectrograph (one of b1,r1,b2,r2) and [eeeeeeee] is the
                zero-padded exposure number. The extension [ext] is "fits" for
                spCFrame files and "fits.gz" for all other file types. Returns None if
                the name is unknown for this band and fiber combination.

        Raises:
            ValueError: one of the inputs is invalid.
        """
        if sequence_number < 0 or sequence_number >= self.num_science_exposures:
            raise ValueError('Invalid sequence number ({0}) must be 0-{1}.'.format(
                sequence_number, self.num_science_exposures - 1))
        if fiber < 1 or fiber > self.num_fibers:
            raise ValueError('Invalid fiber ({}) must be 1-{} for plate {}.'.format(
                fiber, self.num_fibers, self.plate))
        if band not in ('blue', 'red'):
            raise ValueError('Invalid band ({}) must be blue or red.'.format(band))
        if ftype not in ('spCFrame', 'spFrame', 'spFluxcalib', 'spFluxcorr'):
            raise ValueError('Invalid file type ({}) must be one of: '.format(ftype) +
                             'spCFrame, spFrame, spFluxcalib, spFluxcorr.')

        camera = band[0] + str(self.get_spectrograph_index(fiber))
        exposure_info = self.exposures['science'][sequence_number]
        if camera not in exposure_info['SPECS']:
            return None
        exposure_id = exposure_info['EXPID']
        name = '{0}-{1}-{2:08d}.fits'.format(ftype, camera, exposure_id)
        if ftype != 'spCFrame':
            name += '.gz'
        return name


class TraceSet(object):
    """A set of interpolating functions along each trace of a half plate.

    TraceSets use the terminology that x is the pixel coordinate along the nominal
    wavelength direction and y is some quantity to be interpolated as a function
    of x. This implementation is based on the original SDSS IDL code:
    https://trac.sdss3.org/browser/repo/idlutils/trunk/pro/trace/traceset2xy.pro

    Note that red and blue CCDs are handled differently, as described
    `here <https://trac.sdss.org/ticket/880>`_::

        The plan is to switch from 1-phase to 2-phase readout on
        the red CCDs in summer 2010. This will effectively make
        the pixels more uniform, and the flat-fields much better.

        A problem introduced will be that the central two rows will
        each be taller by 1/6 pix. That will flat-field, but there
        will be a discontinuity of 1/3 pix across this point.
        Technically, the PSF will also be different for those pixels,
        and the resulting resolution function.

    Args:
        hdu: fitsio HDU containing the trace set data as a binary table.

    Raises:
        ValueError: Unable to initialize a trace set with this HDU.
    """
    def __init__(self, hdu):
        data = hdu.read()
        if data.shape != (1,):
            raise ValueError('TraceSet HDU has unexpected shape {}.'.format(data.shape))
        data = data[0]
        if data['FUNC'] != 'legendre':
            raise ValueError('TraceSet uses unsupported FUNC "{}"'.format(data['FUNC']))

        # Are we including a jump for the 2-phase readout?
        if 'XJUMPVAL' in hdu.get_colnames():
            self.has_jump = True
            self.xjump_lo = data['XJUMPLO']
            self.xjump_hi = data['XJUMPHI']
            self.xjump_val = data['XJUMPVAL']
        else:
            self.has_jump = False

        self.coefs = data['COEFF']
        if len(self.coefs.shape) != 2:
            raise ValueError('TraceSet coefficients have unexpected shape {}.'.format(
                self.coefs.shape))
        self.ntrace = self.coefs.shape[0]

        xmin, xmax = round(data['XMIN']), round(data['XMAX'])
        self.xmid = 0.5 * (xmin + xmax)
        self.xrange = xmax - xmin
        self.nx = int(self.xrange + 1)

        # Prepare the default xpos used by get_xy()
        one_trace = xmin + np.arange(self.nx)
        self.default_xpos = np.tile(one_trace, (self.ntrace, 1))

    def get_y(self, xpos=None, ignore_jump=False):
        """
        Evaluate the interpolating function for each trace.

        Args:
            xpos(numpy.ndarray): Numpy array of shape (ntrace,nx) with x-pixel
                coordinates along each trace where y(x) should be evaluated. For BOSS,
                ntrace = 500 and for SDSS-I/II (plate < 3510), ntrace = 320. The value
                of ntrace is available as `self.ntrace`.
                If this argument is not set, ``self.default_xpos`` will be used, which
                consists of num_fibers identical traces with x-pixel coordinates at
                each integer pixel value covering the full allowed range.
            ignore_jump(bool): Include a jump when this is set and this is a 2-phase readout.
                There is probably no good reason to set this False, but it is included
                for compatibility with the original IDL code.

        Returns:
            numpy.ndarray: Numpy array ``y`` with shape (ntrace,nx) that matches the input
                ``xpos`` or else the default ``self.default_xpos``.  ``ypos[[i,x]]`` gives
                the value of the interpolated y(x) with x equal to ``xpos[[i,x]]``.
        """
        if xpos is None:
            xpos = self.default_xpos
        y = np.zeros_like(xpos)
        for i in range(self.ntrace):
            x = np.copy(xpos[i])
            if self.has_jump and not ignore_jump:
                t = (x - self.xjump_lo) / (self.xjump_hi - self.xjump_lo)
                below = x < self.xjump_hi
                above = x >= self.xjump_lo
                jump_frac = 1.0 * (~below) + t * (below & above)
                x += jump_frac * self.xjump_val
            xvec = 2 * (x - self.xmid) / self.xrange
            y[i] = numpy.polynomial.legendre.legval(xvec, self.coefs[i])
        return y


class PlateFile(object):
    """A BOSS plate file containing combined exposures for a whole plate.

    This class provides an interface to the :datamodel:`spPlate <PLATE4/spPlate>` data
    product, containing all co-added spectra for a single observation.  To instead
    read individual co-added spectra, use :class:`bossdata.spec.SpecFile`. To access
    individual exposures of a half-plate use :class:`FrameFile`.

    Use :meth:`get_valid_data` to access this plate's data, or the :class:`exposures
    <Exposures>` attribute for a list of exposures used in the coadd.
    The ``num_exposures`` attribute gives the number of science exposures used for this
    target's co-added spectrum (counting a blue+red pair as one exposure). The
    ``plug_map`` attribute records this plate's `plug map
    <http://data.sdss3.org/datamodel/files/PLATELIST_DIR/runs/PLATERUN/plPlugMap.html>`__.

    This class is only intended for reading the BOSS plate file format, so generic
    operations on spectroscopic data (redshifting, resampling, etc) are intentionally not
    included here, but are instead provided in the `speclite
    <http://speclite.readthedocs.org>`__ package.

    Args:
        path(str): Local path of the plate FITS file to use.  This should normally be obtained
            via :meth:`bossdata.path.Finder.get_plate_spec_path` and can be automatically
            mirrored via :meth:`bossdata.remote.Manager.get` or using the :ref:`bossfetch`
            script. The file is opened in read-only mode so you do not need write privileges.
    """
    def __init__(self, path):
        self.hdulist = fitsio.FITS(path, mode=fitsio.READONLY)
        self.header = self.hdulist[0].read_header()
        # Look up the number of fibers.
        self.num_fibers = self.header['NAXIS2']
        # Look up the number of exposures used for this coadd.
        self.exposures = Exposures(self.header)
        self.num_exposures = len(self.exposures.sequence)
        # Calculate the common wavelength grid from header keywords.
        num_pixels = self.header['NAXIS1']
        loglam_min = self.header['COEFF0']
        loglam_step = self.header['COEFF1']
        # We cast the log(lambda) grid to float32 to match how the grid is calculated
        # in IDL and propagated to spec files, even though this entails a (tiny) loss
        # of precision.
        self.loglam = (loglam_min + loglam_step * np.arange(num_pixels)).astype(np.float32)
        # Do not read arrays until we have to.
        self.masks = None
        self.ivar = None
        self.flux = None
        self.wdisp = None
        self.sky = None
        # Read our plug map into an astropy table.
        self.plug_map = astropy.table.Table(self.hdulist[5].read())

    def get_fiber_offsets(self, fiber):
        """Convert fiber numbers to array offsets.

        Args:
            fibers(numpy.ndarray): Numpy array of fiber numbers 1-1000 (or 1-640 for
                plate < 3510).  Fibers do not need to be sorted and repetitions are ok.

        Returns:
            numpy.ndarray: Numpy array of offsets 0-999.

        Raises:
            ValueError: Fiber number is out of the valid range for this plate.
        """
        offset = fiber - 1
        if np.any((offset < 0) | (offset > 999)):
            raise ValueError('Fiber number out of range for this plate.')
        return offset

    def get_pixel_masks(self, fibers):
        """Get the pixel masks for specified fibers.

        The entire mask is returned for each fiber, including any pixels with zero
        inverse variance. Returns the 'and_mask' and ignores the 'or_mask'.

        Args:
            fibers(numpy.ndarray): Numpy array of fiber numbers 1-1000 (or 1-640 for
                plate < 3510).  Fibers do not need to be sorted and repetitions are ok.

        Returns:
            numpy.ndarray: Integer numpy array of shape (nfibers,npixels) where (i,j)
                encodes the mask bits defined in :attr:`bossdata.bits.SPPIXMASK` (see also
                http://www.sdss3.org/dr10/algorithms/bitmask_sppixmask.php) for pixel-j
                of the fiber with index fibers[i].
        """
        offsets = self.get_fiber_offsets(fibers)
        if self.masks is None:
            self.masks = self.hdulist[2].read()
        return self.masks[offsets]

    def get_valid_data(self, fibers, pixel_quality_mask=None,
                       include_wdisp=False, include_sky=False, use_ivar=False,
                       use_loglam=False, fiducial_grid=False):
        """Get the valid for the specified fibers.

        Args:
            fibers(numpy.ndarray): Numpy array of fiber numbers 1-1000 (or 1-640 for
                plate < 3510).  Fibers do not need to be sorted and repetitions are ok.
            pixel_quality_mask(int): An integer value interpreted as a bit pattern using the
                bits defined in :attr:`bossdata.bits.SPPIXMASK` (see also
                http://www.sdss3.org/dr10/algorithms/bitmask_sppixmask.php). Any bits set in
                this mask are considered harmless and the corresponding spectrum pixels are
                assumed to contain valid data. This mask is applied to the AND of the masks
                for each individual exposure. No mask is applied if this value is None.
            include_wdisp: Include a wavelength dispersion column in the returned data.
            include_sky: Include a sky flux column in the returned data.
            use_ivar: Replace ``dflux`` with ``ivar`` (inverse variance) in the returned
                data.
            use_loglam: Replace ``wavelength`` with ``loglam`` (``log10(wavelength)``) in
                the returned data.
            fiducial_grid: Return co-added data using the :attr:`fiducial wavelength grid
                <bossdata.spec.fiducial_loglam>`.  If False, the returned array uses
                the native grid of the SpecFile, which generally trims pixels on both ends
                that have zero inverse variance.  Set this value True to ensure that all
                co-added spectra use aligned wavelength grids when this matters.

        Returns:
            numpy.ma.MaskedArray: Masked array of shape (nfibers,npixels). Pixels with no
                valid data are included but masked. The record for each pixel has at least
                the following named fields: wavelength in Angstroms (or loglam), flux and
                dflux in 1e-17 ergs/s/cm2/Angstrom (or flux and ivar). Wavelength values
                are strictly increasing and dflux is calculated as ivar**-0.5 for pixels
                with valid data. Optional fields are wdisp in constant-log10-lambda pixels
                and sky in 1e-17 ergs/s/cm2/Angstrom. The wavelength (or loglam) field is
                never masked and all other fields are masked when ivar is zero or a
                pipeline flag is set (and not allowed by ``pixel_quality_mask``).
        """
        offsets = self.get_fiber_offsets(fibers)
        num_fibers = len(offsets)

        # Apply the pixel quality mask, if any.
        pixel_bits = self.get_pixel_masks(fibers)
        if pixel_quality_mask is not None:
            clear_allowed = np.bitwise_not(np.uint32(pixel_quality_mask))
            pixel_bits = np.bitwise_and(pixel_bits, clear_allowed)

        # Read arrays from the FITS file if necessary.
        if self.ivar is None:
            self.ivar = self.hdulist[1].read()
        if self.flux is None:
            self.flux = self.hdulist[0].read()
        if include_wdisp and self.wdisp is None:
            self.wdisp = self.hdulist[4].read()
        if include_sky and self.sky is None:
            self.sky = self.hdulist[6].read()

        if fiducial_grid:
            loglam = fiducial_loglam
            first_index = float(get_fiducial_pixel_index(10.0**self.loglam[0]))
            if abs(first_index - round(first_index)) > 0.01:
                raise RuntimeError('Wavelength grid not aligned with fiducial grid.')
            trimmed = slice(first_index, first_index + pixel_bits.shape[1])
        else:
            loglam = self.loglam
            trimmed = slice(None)
        num_pixels = len(loglam)

        # Identify the pixels with valid data.
        ivar = self.ivar[offsets]
        bad_pixels = (pixel_bits != 0) | (ivar <= 0.0)
        good_pixels = ~bad_pixels

        # Create and fill the unmasked structured array of data.
        dtype = [('loglam' if use_loglam else 'wavelength', np.float32),
                 ('flux', np.float32), ('ivar' if use_ivar else 'dflux', np.float32)]
        if include_wdisp:
            dtype.append(('wdisp', np.float32))
        if include_sky:
            dtype.append(('sky', np.float32))
        data = np.zeros((num_fibers, num_pixels), dtype=dtype)
        if use_loglam:
            data['loglam'][:] = loglam
        else:
            data['wavelength'][:] = np.power(10.0, loglam)
        data['flux'][:,trimmed] = self.flux[offsets]
        if use_ivar:
            data['ivar'][:,trimmed][good_pixels] = self.ivar[offsets][good_pixels]
        else:
            data['dflux'][:,trimmed][good_pixels] = 1.0 / np.sqrt(self.ivar[offsets][good_pixels])
        if include_wdisp:
            data['wdisp'][:,trimmed] = self.wdisp[offsets]
        if include_sky:
            data['sky'][:,trimmed] = self.sky[offsets]

        if fiducial_grid:
            mask = np.ones_like(data, dtype=bool)
            mask[:,trimmed] = bad_pixels
        else:
            mask = bad_pixels

        result = numpy.ma.MaskedArray(data, mask=mask)
        # Wavelength values are always valid.
        result['loglam' if use_loglam else 'wavelength'].mask = False

        return result


class FrameFile(object):
    """A BOSS frame file containing a single exposure of one spectrograph (half plate).

    This class supports both types of frame data files: the uncalibrated
    :datamodel:`spFrame <PLATE4/spFrame>` and the calibrated :datamodel:`spCFrame
    <PLATE4/spCFrame>`. Use :meth:`get_valid_data` to access this plate's data and
    the ``plug_map`` attribute to access this plate's `plug map
    <http://data.sdss3.org/datamodel/files/PLATELIST_DIR/runs/PLATERUN/plPlugMap.html>`__.

    BOSS spectrographs read out 500 fibers each. SDSS-I/II spectrographs (plate < 3510)
    read out 320 fibers each.  The ``plate``, ``camera`` and ``exposure_id`` attributes
    provide the basic metadata for this exposure.  The complete HDU0 header is available
    as the ``header`` attribute.

    This class is only intended for reading the BOSS frame file format, so generic
    operations on spectroscopic data (redshifting, resampling, etc) are intentionally not
    included here, but are instead provided in the `speclite
    <http://speclite.readthedocs.org>`__ package.

    Args:
        path(str): Local path of the frame FITS file to use.  This should normally be obtained
            via :meth:`Plan.get_exposure_name` and can be automatically mirrored
            via :meth:`bossdata.remote.Manager.get` or using the :ref:`bossfetch` script. The
            file is opened in read-only mode so you do not need write privileges.
        index(int): Identifies if this is the first (1) or second (2) spectrograph, which
            determines whether it has spectra for fibers 1-500 (1-320) or 501-1000 (321-640).
            You should normally obtain this value using :meth:`Plan.get_spectrograph_index`.
            As of v0.2.7, this argument is optional and will be inferred from the file header
            when not provided, or checked against the file header when provided.
        calibrated(bool): Identifies whether this is a calibrated (spCFrame) or
            un-calibrated (spFrame) frame file. As of v0.2.7, this argument is optional and
            will be inferred from the file header when not provided, or checked against the
            file header when provided.
    """
    def __init__(self, path, index=None, calibrated=None):
        self.hdulist = fitsio.FITS(path, mode=fitsio.READONLY)
        self.header = self.hdulist[0].read_header()
        # Look up the number of fibers.
        self.num_fibers = self.header['NAXIS2']
        # Do not read arrays until we have to.
        self.masks = None
        self.ivar = None
        self.loglam = None
        self.flux = None
        self.wdisp = None
        self.sky = None
        # Read our plug map into an astropy table.
        self.plug_map = astropy.table.Table(self.hdulist[5].read())
        # Extract some metadata from the HDU0 header.
        self.plate = self.header['PLATEID']
        self.exposure_id = self.header['EXPOSURE']
        self.camera = self.header['CAMERAS'].rstrip()
        if self.camera not in ('b1', 'b2', 'r1', 'r2'):
            raise RuntimeError('Found unexpected camera name: {}.'.format(self.camera))
        self.index = int(self.camera[1])
        # Use the HDU3 (wavelength solution) dimensions to detect if this file contains
        # flux-calibrated spectra.
        hdr3 = self.hdulist[3].read_header()
        naxis2 = hdr3['NAXIS2']
        if naxis2 == self.num_fibers:
            self.calibrated = True
        elif naxis2 == 1:
            self.calibrated = False
        else:
            raise RuntimeError('Found unexpected HDU3 NAXIS2 = {}.'.format(naxis2))
        # Check the values of the optional index & calibrated args if they are provided.
        if index is not None and index != self.index:
            raise ValueError('Specified index ({}) does not match the file value ({}).'
                             .format(index, self.index))
        if calibrated is not None and calibrated != self.calibrated:
            raise ValueError('Specified calibrated ({}) does not match the file value ({}).'
                             .format(calibrated, self.calibrated))

    def get_fiber_offsets(self, fiber):
        """Convert fiber numbers to array offsets.

        Args:
            fibers(numpy.ndarray): Numpy array of fiber numbers 1-1000 (or 1-640 for
                plate < 3510).  All fibers must be in the appropriate range 1-500 (1-320)
                or 501-1000 (321-640) for this frame's spectograph. Fibers do not need
                to be sorted and repetitions are ok.

        Returns:
            numpy.ndarray: Numpy array of offsets 0-499 (or 0-319 for plate < 3510).

        Raises:
            ValueError: Fiber number is out of the valid range for this spectrograph.
        """
        offset = fiber - self.num_fibers * (self.index - 1) - 1
        if np.any((offset < 0) | (offset >= self.num_fibers)):
            raise ValueError('Fiber number out of range for this spectrograph.')
        return offset

    def get_pixel_masks(self, fibers):
        """Get the pixel masks for specified fibers.

        The entire mask is returned for each fiber, including any pixels with zero
        inverse variance.

        Args:
            fibers(numpy.ndarray): Numpy array of fiber numbers 1-1000 (or 1-640 for
                plate < 3510).  All fibers must be in the appropriate range 1-500 (1-320)
                or 501-1000 (321-640) for this frame's spectograph. Fibers do not need
                to be sorted and repetitions are ok.

        Returns:
            numpy.ndarray: Integer numpy array of shape (nfibers,npixels) where (i,j)
                encodes the mask bits defined in :attr:`bossdata.bits.SPPIXMASK` (see also
                http://www.sdss3.org/dr10/algorithms/bitmask_sppixmask.php) for pixel-j
                of the fiber with index fibers[i].
        """
        offsets = self.get_fiber_offsets(fibers)
        if self.masks is None:
            self.masks = self.hdulist[2].read()
        return self.masks[offsets]

    def get_valid_data(self, fibers, pixel_quality_mask=None, include_wdisp=False,
        include_sky=False, use_ivar=False, use_loglam=False):
        """Get the valid for the specified fibers.

        Args:
            fibers(numpy.ndarray): Numpy array of fiber numbers 1-1000 (or 1-640 for
                plate < 3510).  All fibers must be in the appropriate range 1-500 (1-320)
                or 501-1000 (321-640) for this frame's spectograph. Fibers do not need
                to be sorted and repetitions are ok.
            pixel_quality_mask(int): An integer value interpreted as a bit pattern using the
                bits defined in :attr:`bossdata.bits.SPPIXMASK` (see also
                http://www.sdss3.org/dr10/algorithms/bitmask_sppixmask.php). Any bits set in
                this mask are considered harmless and the corresponding spectrum pixels are
                assumed to contain valid data.
            include_wdisp: Include a wavelength dispersion column in the returned data.
            include_sky: Include a sky flux column in the returned data.
            use_ivar: Replace ``dflux`` with ``ivar`` (inverse variance) in the returned
                data.
            use_loglam: Replace ``wavelength`` with ``loglam`` (``log10(wavelength)``) in
                the returned data.

        Returns:
            numpy.ma.MaskedArray: Masked array of shape (nfibers,npixels). Pixels with no
                valid data are included but masked. The record for each pixel has at least
                the following named fields: wavelength in Angstroms, flux and dflux in 1e-17
                ergs/s/cm2/Angstrom. Wavelength values are strictly increasing and dflux is
                calculated as ivar**-0.5 for pixels with valid data. Optional fields are
                wdisp in constant-log10-lambda pixels and sky in 1e-17 ergs/s/cm2/Angstrom.
                The wavelength (or loglam) field is never masked and
                all other fields are masked when ivar is zero or a pipeline flag is set (and
                not allowed by ``pixel_quality_mask``).
        """
        offsets = self.get_fiber_offsets(fibers)
        num_fibers = len(offsets)

        # Apply the pixel quality mask, if any.
        pixel_bits = self.get_pixel_masks(fibers)
        if pixel_quality_mask is not None:
            clear_allowed = np.bitwise_not(np.uint32(pixel_quality_mask))
            pixel_bits = np.bitwise_and(pixel_bits, clear_allowed)

        # Read arrays from the FITS file if necessary.
        if self.ivar is None:
            self.ivar = self.hdulist[1].read()
        if self.loglam is None:
            if self.calibrated:
                self.loglam = self.hdulist[3].read()
            else:
                # Expand the traceset solution. Wavelengths are stored as log10(lambda).
                self.wavelength_trace_set = TraceSet(self.hdulist[3])
                self.loglam = self.wavelength_trace_set.get_y()
                if self.loglam.shape != self.ivar.shape:
                    raise RuntimeError('HDU3 traceset has unexpected shape: {}.'.format(
                        self.loglam.shape))
        if self.flux is None:
            self.flux = self.hdulist[0].read()
        if include_wdisp and self.wdisp is None:
            if self.calibrated:
                self.wdisp = self.hdulist[4].read()
            else:
                # Expand the traceset solution. Dispersions are in units of 1e4*d(lambda)/dx.
                xpos = self.wavelength_trace_set.default_xpos
                loglam_hi = self.wavelength_trace_set.get_y(xpos + 0.5)
                loglam_lo = self.wavelength_trace_set.get_y(xpos - 0.5)
                dlambda_by_dx = np.abs(loglam_hi - loglam_lo)
                dispersion_trace_set = TraceSet(self.hdulist[4])
                self.wdisp = 1e4 * dlambda_by_dx * dispersion_trace_set.get_y()
                if self.wdisp.shape != self.ivar.shape:
                    raise RuntimeError('HDU4 traceset has unexpected shape: {}.'.format(
                        self.wdisp.shape))
        if include_sky and self.sky is None:
            self.sky = self.hdulist[6].read()
        num_pixels = self.flux.shape[1]

        # Identify the pixels with valid data.
        ivar = self.ivar[offsets]
        bad_pixels = (pixel_bits != 0) | (ivar <= 0.0)
        good_pixels = ~bad_pixels

        # Create and fill the unmasked structured array of data.
        dtype = [('loglam' if use_loglam else 'wavelength', np.float32),
                 ('flux', np.float32), ('ivar' if use_ivar else 'dflux', np.float32)]
        if include_wdisp:
            dtype.append(('wdisp', np.float32))
        if include_sky:
            dtype.append(('sky', np.float32))
        data = np.empty((num_fibers, num_pixels), dtype=dtype)
        if use_loglam:
            data['loglam'][:] = self.loglam[offsets]
        else:
            data['wavelength'][:] = np.power(10.0, self.loglam[offsets])
        data['flux'][:] = self.flux[offsets]
        if use_ivar:
            data['ivar'][:][good_pixels] = self.ivar[offsets][good_pixels]
        else:
            data['dflux'][:][good_pixels] = 1.0 / np.sqrt(self.ivar[offsets][good_pixels])
        if include_wdisp:
            data['wdisp'][:] = self.wdisp[offsets]
        if include_sky:
            data['sky'][:] = self.sky[offsets]

        mask = bad_pixels
        result = numpy.ma.MaskedArray(data, mask=mask)
        # Wavelength values are always valid.
        result['loglam' if use_loglam else 'wavelength'].mask = False

        return result

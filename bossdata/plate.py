# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Access BOSS plate data products.
"""

from __future__ import division, print_function

import os.path

import numpy as np
import numpy.ma
import numpy.polynomial.legendre

import fitsio


class Plan(object):
    """The plan file for configuring the BOSS pipeline to combine exposures of a single plate.

    The datamodel for plan files is at
    http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/PLATE4/spPlan.html.
    Combined plan files are small text files that list the per-spectrograph (b1,b2,r1,r2)
    spFrame files used for a single coadd.

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
                    self.plate = tokens[1]
                elif self.plate != tokens[1]:
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

    def get_exposure_name(self, sequence_number, camera, fiber, calibrated=True):
        """Get the name of a science exposure from this plan.

        Use the exposure name to locate the calibrated spCFrame-{EXPNAME}.fits and
        uncalibrated spFrame-{EXPNAME}.fits.gz files for individual exposures.

        Args:
            sequence_number(int): Science exposure sequence number, counting from zero.
                Must be less than our num_science_exposures attribute.
            fiber(int): Fiber number to identify which spectrograph to use, which must
                be in the range 1-1000.
            camera(str): Must be 'blue' or 'red'.
            calibrated(bool): Returns the name of the calibrated (spCFrame) file rather
                than the un-calibrated (spFrame) file.

        Returns:
            str: Exposure name of the form [prefix]-[cc]-[eeeeeeee].[ext] where [cc]
                identifies the spectrograph (one of b1,r1,b2,r2) and [eeeeeeee] is the
                zero-padded exposure number. For calibrated exposures, [prefix] is
                "spCFrame" and [ext] is "fits".  For un-calibrated exposures, [prefix]
                is "spFrame" and [ext] is "fits.gz". Returns None if the name is
                unknown for this camera and fiber combination.

        Raises:
            ValueError: one of the inputs is invalid.
        """
        if sequence_number < 0 or sequence_number >= self.num_science_exposures:
            raise ValueError('Invalid sequence number ({0}) must be 0-{1}.'.format(
                sequence_number, self.num_science_exposures))
        if fiber < 1 or fiber > 1000:
            raise ValueError('Invalid fiber ({}) must be 1-1000.'.format(fiber))
        if camera not in ('blue', 'red'):
            raise ValueError('Invalid camera ({}) must be blue or red.'.format(camera))

        if fiber <= 500:
            index = '1'
        else:
            index = '2'
        spectrograph = camera[0] + index
        exposure_info = self.exposures['science'][sequence_number]
        if spectrograph not in exposure_info['SPECS']:
            return None
        exposure_id = exposure_info['EXPID']
        if calibrated:
            prefix, ext = 'spCFrame', 'fits'
        else:
            prefix, ext = 'spFrame', 'fits.gz'

        return '{0}-{1}-{2:08d}.{3}'.format(prefix, spectrograph, exposure_id, ext)


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
            xpos(numpy.ndarray): Numpy array of shape (500,nx) with x-pixel coordinates
                along each trace where y(x) should be evaluated. If this
                argument is not set, ``self.default_xpos`` will be used, which consists
                of 500 identical traces with x-pixel coordinates at each integer pixel
                value covering the full allowed range (nominally 0-4111 for blue, 0-4127
                for red).
            ignore_jump(bool): Include a jump when this is set and this is a 2-phase readout.
                There is probably no good reason to set this False, but it is included
                for compatibility with the original IDL code.

        Returns:
            numpy.ndarray: Numpy array ``y`` with shape (500,nx) that matches the input
                ``xpos`` or else the default ``self.default_xpos``.  ``ypos[[i,x]]`` gives
                the value of the interpolated y(x) with x equal to ``xpos[[i,x]]``.
        """
        if xpos is None:
            xpos = self.default_xpos
        y = np.zeros_like(xpos)
        for i in range(self.ntrace):
            x = np.copy(xpos[i])
            if self.has_jump and not ignore_jump:
                t = (x - self.xjump_lo)/(self.xjump_hi - self.xjump_lo)
                below = x < self.xjump_hi
                above = x >= self.xjump_lo
                jump_frac = 1.0 * (~below) + t * (below & above)
                x += jump_frac * self.xjump_val
            xvec = 2 * (x - self.xmid) / self.xrange
            y[i] = numpy.polynomial.legendre.legval(xvec, self.coefs[i])
        return y

class FrameFile(object):
    """A BOSS frame file containing a single exposure of one spectrograph (500 fibers).

    This class supports both types of frame data files: the uncalibrated spFrame and
    the calibrated spCFrame. The corresponding data models are documented at:

    http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/PLATE4/spFrame.html
    http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/PLATE4/spCFrame.html

    Args:
        path(str): Local path of the frame FITS file to use.  This should normally be obtained
            via :meth:`Plan.get_exposure_name` and can be automatically mirrored
            via :meth:`bossdata.remote.Manager.get` or using the :ref:`bossfetch` script. The
            file is opened in read-only mode so you do not need write privileges.
        index(int): Identifies if this is the first (1) or second (2) spectrograph, which
            determines whether it has spectra for fibers 1-500 or 501-1000.
        calibrated(bool): Identifies whether this is a calibrated (spCFrame) or
            un-calibrated (spFrame) frame file.
    """
    def __init__(self, path, index, calibrated):
        if index not in (1, 2):
            raise ValueError('Invalid index ({}) should be 1 or 2.'.format(index))
        self.index = index
        self.calibrated = calibrated
        self.hdulist = fitsio.FITS(path, mode=fitsio.READONLY)
        self.masks = None
        self.ivar = None
        self.loglam = None
        self.flux = None
        self.wdisp = None
        self.sky = None

    def get_fiber_offsets(self, fiber):
        """Convert fiber numbers to array offsets.

        Args:
            fibers(numpy.ndarray): Numpy array of fiber numbers 1-1000.  All fibers must
                be in the appropriate range 1-500 or 501-1000 for this frame's spectograph.
                Fibers do not need to be sorted and repetitions are ok.

        Returns:
            numpy.ndarray: Numpy array of offsets 0-499.

        Raises:
            ValueError: Fiber number is out of the valid range for this spectrograph.
        """
        offset = fiber - 500*(self.index-1) - 1
        if np.any((offset < 0) | (offset > 499)):
            raise ValueError('Fiber number out of range for this spectrograph.')
        return offset

    def get_pixel_masks(self, fibers):
        """Get the pixel masks for specified fibers.

        The entire mask is returned for each fiber, including any pixels with zero
        inverse variance.

        Args:
            fibers(numpy.ndarray): Numpy array of fiber numbers 1-1000.  All fibers must
                be in the appropriate range 1-500 or 501-1000 for this frame's spectograph.
                Fibers do not need to be sorted and repetitions are ok.

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
                       include_wdisp=False, include_sky=False):
        """Get the valid for a specified exposure or the combined coadd.

        Args:
            fibers(numpy.ndarray): Numpy array of fiber numbers 1-1000.  All fibers must
                be in the appropriate range 1-500 or 501-1000 for this frame's spectograph.
                Fibers do not need to be sorted and repetitions are ok.
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
            numpy.ma.MaskedArray: Masked array of shape (nfibers,npixels). Pixels with no
                valid data are included but masked. The record for each pixel has at least
                the following named fields: wavelength in Angstroms, flux and dflux in 1e-17
                ergs/s/cm2/Angstrom. Wavelength values are strictly increasing and dflux is
                calculated as ivar**-0.5 for pixels with valid data. Optional fields are
                wdisp in constant-log10-lambda pixels and sky in 1e-17 ergs/s/cm2/Angstrom.
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
        dtype = [('wavelength', np.float32), ('flux', np.float32), ('dflux', np.float32)]
        if include_wdisp:
            dtype.append(('wdisp', np.float32))
        if include_sky:
            dtype.append(('sky', np.float32))
        data = np.empty((num_fibers, num_pixels), dtype=dtype)
        data['wavelength'][:] = np.power(10.0, self.loglam[offsets])
        data['flux'][:] = self.flux[offsets]
        data['dflux'][:][good_pixels] = 1.0 / np.sqrt(self.ivar[offsets][good_pixels])
        if include_wdisp:
            data['wdisp'][:] = self.wdisp[offsets]
        if include_sky:
            data['sky'][:] = self.sky[offsets]

        return numpy.ma.MaskedArray(data, mask=bad_pixels)

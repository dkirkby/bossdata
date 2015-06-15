# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Access BOSS plate data products.
"""

from __future__ import division, print_function

import os.path

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
        with open(path,'r') as f:
            for line in f:
                if not line.startswith('SPEXP'):
                    continue
                tokens = line.split()
                # token[1] is the plate number and must be identical for all lines.
                if not self.plate:
                    self.plate = tokens[1]
                elif self.plate != tokens[1]:
                    raise RuntimeError('Internal error: unexpected plate {0} in {1}.'.format(
                        tokens[1],path))
                # Validate the exposure names, which should all have the form
                # [prefix]-[cc]-[eeeeeeee].fits
                exposure_id = set()
                specs = set()
                for name in tokens[7:11]:
                    if name == 'UNKNOWN':
                        continue
                    name,ext = os.path.splitext(name)
                    if ext != '.fits':
                        raise RuntimeError('Unexpected extension {} in {}.'.format(ext,path))
                    fields = name.split('-')
                    if len(fields) != 3:
                        raise RuntimeError('Unexpected exposure name {} in {}.'.format(name,path))
                    if fields[0] not in ('sdR','spFrame'):
                        raise RuntimeError('Unexpected prefix {} in {}.'.format(fields[0],path))
                    if fields[1] not in ('r1','b1','r2','b2'):
                        raise RuntimeError('Unexpected camera {} in {}.'.format(fields[1],path))
                    exposure_id.add(int(fields[2]))
                    specs.add(fields[1])
                if len(exposure_id) != 1:
                    raise RuntimeError('Multiple exposure IDs: {}.'.format(exposure_id))
                # Build an exposure record to save.
                exposure = dict(
                    MJD=tokens[2], EXPTIME=float(tokens[5]), EXPID=exposure_id.pop(), SPECS=specs)
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
            calibrated(bool): Returns the name of the flux-calibrated exposure.

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
        if camera not in ('blue','red'):
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
            prefix, ext = 'spCFrame','fits'
        else:
            prefix, ext = 'spFrame','fits.gz'

        return '{0}-{1}-{2:08d}.{3}'.format(prefix,spectrograph,exposure_id,ext)

class Frame(object):
    """A single-exposure frame of one BOSS spectrograph.

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
    """
    def __init__(self, path, index):
        if index not in (1,2):
            raise ValueError('Invalid index ({}) should be 1 or 2.'.format(index))
        self.index = index
        self.hdulist = fitsio.FITS(path, mode=fitsio.READONLY)
        self.masks = None

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
        offset = fiber - 500*(index+1) - 1
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

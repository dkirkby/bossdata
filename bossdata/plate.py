# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Access BOSS plate data products.
"""

from __future__ import division, print_function

import os.path

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
                for name in tokens[7:11]:
                    name,ext = os.path.splitext(name)
                    if ext != '.fits':
                        raise RuntimeError('Unexpected extension {}.'.format(ext))
                    fields = name.split('-')
                    if len(fields) != 3:
                        raise RuntimeError('Unexpected exposure name {}.'.format(name))
                    if fields[0] not in ('sdR','spFrame'):
                        raise RuntimeError('Unexpected prefix {}.'.format(fields[0]))
                    if fields[1] not in ('r1','b1','r2','b2'):
                        raise RuntimeError('Unexpected camera {}.'.format(fields[1]))
                    exposure_id.add(int(fields[2]))
                if len(exposure_id) != 1:
                    raise RuntimeError('Multiple exposure IDs: {}.'.format(exposure_id))
                # Build an exposure record to save.
                exposure = dict(
                    MJD=tokens[2], EXPTIME=float(tokens[5]), EXPID=exposure_id.pop())
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
                is "spFrame" and [ext] is "fits.gz".

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
        exposure_id = self.exposures['science'][sequence_number]['EXPID']
        if calibrated:
            prefix, ext = 'spCFrame','fits'
        else:
            prefix, ext = 'spFrame','fits.gz'

        return '{0}-{1}-{2:08d}.{3}'.format(prefix,spectrograph,exposure_id,ext)

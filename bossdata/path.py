# -*- coding: utf-8 -*-

"""Generate paths to BOSS data files.
"""

from __future__ import division,print_function

import os
import os.path

class Finder(object):
    """Initialize a path finder object.

    Args:
        sas_root(str): Location of the SAS root path to use, e.g., /sas/dr12. Will use the
            value of the BOSS_SAS_ROOT environment variable if this is not set.
        redux_version(str): String tag specifying the BOSS spectro reduction version to use,
            e.g., v5_7_0. Will use the value of the BOSS_REDUX_VERSION environment variable
            if this is not set.

    Raises:
        RuntimeError: No SAS root or redux version specified on the command line or via
            environment variables.
    """
    def __init__(self,sas_root=None,redux_version=None):
        if sas_root is None:
            sas_root = os.getenv('BOSS_SAS_ROOT')
        if sas_root is None:
            raise RuntimeError('No SAS root specified: try setting $BOSS_SAS_ROOT.')
        if redux_version is None:
            redux_version = os.getenv('BOSS_REDUX_VERSION',None)
        if redux_version is None:
            raise RuntimeError('No redux version specifed: try setting $BOSS_REDUX_VERSION.')

        self.redux_version = redux_version
        self.redux_base = os.path.join(sas_root,'boss','spectro','redux',redux_version)

    def get_plate_path(self,plate):
        """Get the path to the specified plate.

        The returned path contains files that include all targets on the plate. Use the
        :meth:`get_spec_path` method for the path of a single spectrum file.

        This method only performs minimal checks that the requested plate number is valid.

        Args:
            plate(int): Plate number, which must be positive.

        Returns:
            str: Full path to files for the specified plate.

        Raises:
            ValueError: Invalid plate number must be > 0.
        """
        if plate < 0:
            raise ValueError('Invalid plate number ({}) must be > 0.'.format(plate))

        return os.path.join(self.redux_base,str(plate))

    def get_sp_all_path(self):
        """Get the location of the metadata summary file.
        """
        name = 'spAll-{}.dat.gz'.format(self.redux_version)
        return os.path.join(self.redux_base,name)

    def get_spec_path(self,plate,mjd,fiber):
        """Get the location of the spectrum file for the specified observation.

        The DR12 data model for the returned files is at
        http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html.
        Each file is approximately 1.5Mb in size.

        Use the :meth:`get_plate_path` method for the path to files that include all targets
        on a plate.

        This method only performs minimal checks that the requested plate-mjd-fiber are valid.

        Args:
            plate(int): Plate number, which must be positive.
            mjd(int): Modified Julian date of the observation, which must be > 3500.
            fiber(int): Fiber number of the target on this plate, which must be in the range 1-1000.

        Returns:
            str: Full path to the spectrum file for the specified observation.

        Raises:
            ValueError: Invalid plate number must be > 0.
        """
        if plate < 0:
            raise ValueError('Invalid plate number ({}) must be > 0.'.format(plate))
        if mjd <= 3500:
            raise ValueError('Invalid mjd ({}) must be >= 3500.'.format(mjd))
        if fiber < 1 or fiber > 1000:
            raise ValueError('Invalid fiber ({}) must be 1-1000.'.format(fiber))

        name = 'spec-{plate:4d}-{mjd:5d}-{fiber:04d}.fits'.format(plate=plate,mjd=mjd,fiber=fiber)
        return os.path.join(self.redux_base,'spectra',str(plate),name)

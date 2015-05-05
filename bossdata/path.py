# -*- coding: utf-8 -*-

"""Brief description.

A more detailed description goes here...
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

        self.redux_base = os.path.join(sas_root,'boss','spectro','redux',redux_version)

    def get_plate_path(self,plate):
        """Get the path to the specified plate.

        This method only performs minimal plates that the requested plate number is valid.

        Args:
            plate(int): Plate number, which must be positive.

        Returns:
            str: Full path to files for the specified plate.

        Raises:
            ValueError: Invalid plate number must be > 0.
        """
        if plate < 0:
            raise ValueError('Invalid Invalid plate number ({}) must be > 0.'.format(plate))

        return os.path.join(self.redux_base,str(plate))

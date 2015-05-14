# -*- coding: utf-8 -*-

""" Access spectroscopic data for a single BOSS target.
"""

from __future__ import division,print_function

import fitsio

class SpecFile(object):
    """ A BOSS spec file containing summary data for a single target.

    The spec file data model is described at
    http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html.
    This class supports the full version described in the data model as well as a "lite"
    version that does not contain the per-exposure HDUs with indices >= 4.

    This class is designed only to interface with the BOSS spec file format and generic
    operations on spectroscopic data are deliberately not included here.

    Args:
        path(str): Local path of the spec FITS file to use.  This should normally be obtained
            via :meth:`bossdata.path.get_spec_path` and can be automatically mirrored via
            :meth:`bossdata.remote.get` or using the :ref:`bossfetch` script.
    """
    def __init__(self,path):
        pass

# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

""" Access spectroscopic data for a single BOSS target.
"""

from __future__ import division,print_function

import numpy as np

import fitsio

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
        # Check for non-lite HDUs
        self.lite = (len(self.hdulist) == 4)

    def get_spectrum_hdu(self,exposure_index = None):
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

    def get_wavelength(self,exposure_index = None):
        """Get the wavelength grid of this spectrum.

        Use :meth:`get_flux` and :meth:`get_ivar` to get the flux measurements associated
        with this grid.

        Args:
            exposure_index(int): Individual exposure to use. Uses the co-added spectrum when
                this is None.

        Returns:
            numpy.ndarray: Array of increasing wavelength values in Angstroms.
        """
        log10_wavelength = self.get_spectrum_hdu(exposure_index)['loglam'][:]
        return np.power(10.0,log10_wavelength)

    def get_flux(self,exposure_index = None):
        """Get the flux values for this spectrum.

        Args:
            exposure_index(int): Individual exposure to use. Uses the co-added spectrum when
                this is None.

        Returns:
            numpy.ndarray: Array of flux values in 1e-17 ergs/s/cm2/Angstrom, corresponding to
                the wavelengths returned by :meth:`get_wavelength`.
        """
        return self.get_spectrum_hdu(exposure_index)['flux'][:]

    def get_ivar(self,exposure_index = None):
        """Get the flux inverse variances for this spectrum.

        When ivar is non-zero, the corresponding flux standard deviation is ivar**-0.5.
        The special value ivar=0 is used to indicate that no flux information is available.

        Args:
            exposure_index(int): Individual exposure to use. Uses the co-added spectrum when
                this is None.

        Returns:
            numpy.ndarray: Array of flux inverse variances in 1e-17 ergs/s/cm2/Angstrom,
                corresponding to the wavelengths returned by :meth:`get_wavelength`.
        """
        return self.get_spectrum_hdu(exposure_index)['ivar'][:]

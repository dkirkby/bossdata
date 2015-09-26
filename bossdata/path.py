# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Generate paths to BOSS data files.

The path module provides convenience methods for building the paths of frequently used
data files.  Most scripts will create a single :class:`Finder` object using the default
constructor for this purpose::

    import bossdata.path
    finder = bossdata.path.Finder()

This finder object is normally configured by the `$BOSS_SAS_PATH` and `$BOSS_REDUX_VERSION`
environment variables and no other modules uses these variables, except through a
a :class:`Finder` object.  These parameters can also be set by :class:`Finder` constructor
arguments. When neither the environment variables nor the constructor arguments are set,
defaults appropriate for the most recent public data release (DR12) are used.

:class:`Finder` objects never interact with any local or
remote filesystems: use the :mod:`bossdata.remote` module to download data files
and access them locally. See :doc:`/usage` for recommendations
on using the :mod:`bossdata.path` and :mod:`bossdata.remote` modules together.
"""

from __future__ import division, print_function

import os
import posixpath

from bossdata.plate import get_num_fibers


class Finder(object):
    """Initialize a path finder object.

    When the constructor is called with no arguments, it will raise a ValueError if either
    BOSS_SAS_PATH or BOSS_REDUX_VERSION is not set.

    Args:
        sas_path(str): Location of the SAS root path to use, e.g., /sas/dr12. Will use the
            value of the BOSS_SAS_PATH environment variable if this is not set.
        redux_version(str): String tag specifying the BOSS spectro reduction version to use,
            e.g., v5_7_0. Will use the value of the BOSS_REDUX_VERSION environment variable
            if this is not set.

    Raises:
        ValueError: No SAS root or redux version specified on the command line or via
            environment variables.
    """
    def __init__(self, sas_path=None, redux_version=None, verbose=True):
        # Constructor args take precendence over constructor args.
        self.sas_path = os.getenv('BOSS_SAS_PATH') if sas_path is None else sas_path
        if self.sas_path is None:
            self.sas_path = Finder.default_sas_path
            if verbose:
                print('Using the default "{}" since $BOSS_SAS_PATH is not set.'.format(
                    self.sas_path))

        self.redux_version = (os.getenv('BOSS_REDUX_VERSION') if redux_version is None
                              else redux_version)
        if self.redux_version is None:
            self.redux_version = Finder.default_redux_version
            if verbose:
                print('Using the default "{}" since $BOSS_REDUX_VERSION is not set.'.format(
                    self.redux_version))

        self.redux_base = posixpath.join(self.sas_path, 'spectro', 'redux', self.redux_version)

    default_sas_path = '/sas/dr12/boss'
    """Default to use when $BOSS_SAS_PATH is not set.

    See :doc:`/scripts` and :doc:`/usage` for details.
    """

    default_redux_version = 'v5_7_0'
    """Default to use when $BOSS_REDUX_VERSION is not set.

    See :doc:`/scripts` and :doc:`/usage` for details.
    """

    def get_plate_path(self, plate, filename=None):
        """Get the path to the specified plate directory or file.

        The returned path contains files that include all targets on the plate. Use the
        :meth:`get_spec_path` method for the path of a single spectrum file.

        This method only performs minimal checks that the requested plate number is valid.

        Args:
            plate(int): Plate number, which must be positive.
            filename(str): Name of a file within the plate directory to append to the
                returned path.

        Returns:
            str: Full path to the specified plate directory or file within this directory.

        Raises:
            ValueError: Invalid plate number must be > 0.
        """
        if plate < 0:
            raise ValueError('Invalid plate number ({}) must be > 0.'.format(plate))

        path = posixpath.join(self.redux_base, '{:04d}'.format(plate))
        if filename:
            path = posixpath.join(path, filename)
        return path

    def get_plate_plan_path(self, plate, mjd, combined=True):
        """Get the path to the specified plate plan file.

        A combined plan may span several nearby MJDs, in which case the last MJD is
        the one used to identify the plan.

        Args:
            plate(int): Plate number, which must be positive.
            mjd(int): Modified Julian date of the observation, which must be > 45000.
            combined(bool): Specifies the combined plan, which spans all MJDs
                associated with a coadd, but does not include calibration frames
                (arcs,flats) for a specific MJD.

        Returns:
            str: Full path to the requested plan file.

        Raises:
            ValueError: Invalid plate or mjd inputs.
        """
        if mjd <= 45000:
            raise ValueError('Invalid mjd ({}) must be >= 45000.'.format(mjd))

        if combined:
            filename = 'spPlancomb-{plate:04d}-{mjd:5d}.par'.format(plate=plate, mjd=mjd)
        else:
            filename = 'spPlan2d-{plate:04d}-{mjd:5d}.par'.format(plate=plate, mjd=mjd)
        return posixpath.join(self.get_plate_path(plate), filename)

    def get_plate_spec_path(self, plate, mjd):
        """Get the path to the file containing combined spectra for a whole plate.

        Combined spectra for all exposures of a plate are packaged in
        :datamodel:`spPlate files <spPlate>`. As of DR12, these files are about 110Mb
        for 1000 spectra.

        Args:
            plate(int): Plate number, which must be positive.
            mjd(int): Modified Julian date of the observation, which must be > 45000.

        Returns:
            str: Full path to the requested plan file.

        Raises:
            ValueError: Invalid plate or mjd inputs.
        """
        if mjd <= 45000:
            raise ValueError('Invalid mjd ({}) must be >= 45000.'.format(mjd))

        filename = 'spPlate-{plate:04d}-{mjd:5d}.fits'.format(plate=plate, mjd=mjd)
        return posixpath.join(self.get_plate_path(plate), filename)

    def get_sp_all_path(self, lite=True):
        """Get the location of the metadata summary file.

        The :datamodel:`spAll <spAll>` file provides extensive metadata for all survey
        targets as a FITS file.  There is also a smaller "lite" version containing a
        subset of this metadata in compressed text format. As of DR12, the full file
        size is about 10Gb and the lite file is about 115Mb.

        Args:
            lite(bool): Specifies the "lite" version which contains all rows but only the
                most commonly used subset of columns. The lite version is a compressed (.gz)
                text data file, while the full version is a FITS file.
        """
        if lite:
            name = 'spAll-{}.dat.gz'.format(self.redux_version)
        else:
            name = 'spAll-{}.fits'.format(self.redux_version)
        return posixpath.join(self.redux_base, name)

    def get_platelist_path(self):
        """Get the location of the platelist summary file.

        The :datamodel:`platelist <platelist>` contains one row per observation (PLATE-MJD),
        unlike most other sources of metadata which contain one row per target
        (PLATE-MJD-FIBER).
        """
        name = 'platelist.fits'
        return posixpath.join(self.redux_base, name)

    default_quasar_catalog_name = 'DR12Q'
    """Default quasar catalog name.

    For more info about the BOSS quasar catalog, see
    http://www.sdss.org/dr12/algorithms/boss-dr12-quasar-catalog/
    """

    def get_quasar_catalog_path(self, catalog_name=None):
        """Get the location of the quasar catalog file.

        The `quasar catalog
        <http://data.sdss3.org/datamodel/files/BOSS_QSO/DR12Q/DR12Q.html>`__
        is documented at
        http://www.sdss.org/dr12/algorithms/boss-dr12-quasar-catalog/.
        As of DR12, the file size is about 513Mb.

        Args:
            catalog_name(str): BOSS quasar catalog name. Will use the
                :meth:`get_default_quasar_catalog_name` method if this is not set.
        """
        if catalog_name is None:
            catalog_name = self.default_quasar_catalog_name
        filename = '{}.fits'.format(catalog_name)
        return posixpath.join(self.sas_path, 'qso', catalog_name, filename)

    def get_spec_path(self, plate, mjd, fiber, lite=True):
        """Get the location of the spectrum file for the specified observation.

        The DR12 data model for the returned files is at
        http://dr12.sdss3.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html
        but only HDUs 0-3 are included in the (default) lite format.
        Each lite (full) file is approximately 0.2Mb (1.7Mb) in size.

        Use the :meth:`get_plate_path` method for the path to files that include all targets
        on a plate.

        This method only performs minimal checks that the requested plate-mjd-fiber are valid.

        Args:
            plate(int): Plate number, which must be positive.
            mjd(int): Modified Julian date of the observation, which must be > 45000.
            fiber(int): Fiber number of the target on this plate, which must be in the
                range 1-1000 (or 1-640 for plate < 3510).
            lite(bool): Specifies the "lite" version which contains only HDUs 0-3, so no
                per-exposure data is included.

        Returns:
            str: Full path to the spectrum file for the specified observation.

        Raises:
            ValueError: Invalid plate, mjd or fiber inputs.
        """
        if plate < 0:
            raise ValueError('Invalid plate number ({}) must be > 0.'.format(plate))
        if mjd <= 45000:
            raise ValueError('Invalid mjd ({}) must be >= 45000.'.format(mjd))
        if fiber < 1 or fiber > get_num_fibers(plate):
            raise ValueError('Invalid fiber ({}) must be 1-{} for plate {}.'.format(
                fiber, get_num_fibers(plate), plate))

        path = 'spectra'
        plate_label = '{:04d}'.format(plate)
        if lite:
            path = posixpath.join(path, 'lite')
        name = 'spec-{plate}-{mjd:5d}-{fiber:04d}.fits'.format(
            plate=plate_label, mjd=mjd, fiber=fiber)
        return posixpath.join(self.redux_base, path, plate_label, name)

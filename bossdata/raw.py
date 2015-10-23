# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Work with raw spectroscopic science and calibration images.
"""

import os.path

import numpy as np

import fitsio

import pydl.pydlutils.yanny


class RawImageFile(object):
    """Wrapper for a raw science or calibration sdR image format.

    Each camera (b1,b2,r1,r2) has its own raw image file for a given exposure
    ID. See the `sdR datamodel
    <http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_DATA/MJD/sdR.html>`__
    for details.  The FITS file is opened, read, and closed by the constructor.

    The raw files used in a co-add can be located and opened using
    :meth:`bossdata.spec.SpecFile.get_raw_image`.

    Args:
        path(str): Local path of the FITS image file to use.

    Attributes:
        header(dict): Dictionary of FITS header keyword-value pairs.
        plate(int): Plate number used to record this raw data.
        exposure_id(int): Exposure ID for this raw data.
        camera(str): One of b1, b2, r1, r2.
        flavor(str): One of science, arc, flat.
        data(numpy.ndarray): 2D array of raw pixel values.
    """
    def __init__(self, path):
        hdulist = fitsio.FITS(path, mode=fitsio.READONLY)
        self.header = hdulist[0].read_header()
        self.plate = self.header['PLATEID']
        self.exposure_id = self.header['EXPOSURE']
        self.camera = self.header['CAMERAS'].rstrip()
        self.flavor = self.header['FLAVOR'].rstrip()
        self.data = hdulist[0].read()
        hdulist.close()

    def get_amplifier_region(self, amplifier_index, region_type):
        """Get overscan and data regions for one amplifier.

        Region definitions are taken from the `sdR datamodel
        <http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_DATA/MJD/sdR.html>`__.
        Each amplifier reads out one quadrant of the sensor.  Amplfiers 0-3 have
        outside corner pixel indices [0, 0], [0, -1], [-1, 0], [-1, -1],
        respectively. Results are only valid for MJD >= 55113 (9-Oct-2009).

        Args:
            amplifier_index(int): Amplifier's are indexed 0-3.
            region_type(str): One of data, overscan.

        Returns:
            tuple: Tuple (rows, cols) of pixel slices that define the requested
                region.  The return value can be used to create a view of our
                raw data as ``self.data[rows, cols]``.
        """
        if self.header['MJD'] < 55113:
            raise ValueError('Amplifier regions not available for MJD < 55113.')
        if amplifier_index < 0 or amplifier_index > 3:
            raise ValueError('Invalid amplifier_index {}. Should be 0-3.'.format(
                amplifier_index))
        if region_type not in ('data', 'overscan'):
            raise ValueError('Invalid region type "{}".  Should be data or overscan.')

        amp_row, amp_col = (amplifier_index // 2), (amplifier_index % 2)
        # Add one to end index relative to the slices in the data model since
        # those use IDL notation, where the end index is included in the result.
        if self.camera[0] == 'b':
            rows = slice(56, 2112) if amp_row == 0 else slice(2112, 4168)
            if region_type == 'overscan':
                cols = slice(10, 68) if amp_col == 0 else slice(4284, 4341)
            else:
                cols = slice(128, 2176) if amp_col == 0 else slice(2176, 4224)
        else:
            rows = slice(48, 2112) if amp_row == 0 else slice(2112, 4176)
            if region_type == 'overscan':
                cols = slice(10, 101) if amp_col == 0 else slice(4250, 4341)
            else:
                cols = slice(119, 2176) if amp_col == 0 else slice(2176, 4233)

        return rows, cols

    def get_amplifier_bias(self, amplifier_index, percentile_cut=1):
        """Estimate the amplifier bias.

        Estimate the bias in one amplifier (quadrant) using the truncated mean
        of pixel values in its :meth:`overscan region <get_amplifier_region>`.

        Args:
            amplifier_index(int): Amplifier's are indexed 0-3.
            percentile_cut(float): Percentage of outliers to ignore on both
                sides of the distribution.

        Returns:
            float: Estimated bias.
        """
        if percentile_cut < 0 or percentile_cut >= 50:
            raise ValueError(
                'Invalid percentile cut {}. Expected 0-50.'.format(percentile_cut))
        overscan_data = self.data[
            self.get_amplifier_region(amplifier_index, 'overscan')].flatten()
        lo, hi = np.percentile(overscan_data, (percentile_cut, 100 - percentile_cut))
        in_range = (overscan_data > lo) & (overscan_data < hi)
        return np.mean(overscan_data[in_range])

    def get_data(self, bias_subtracted=True, percentile_cut=1, bias_point=100):
        """Get the data region of this raw image.

        The data region is the union of the four :meth:`amplifier data regions
        <get_amplifier_region>`.

        Args:
            bias_subtracted(bool): Subtract bias from each amplifier quadrant,
                estimated using :meth:`get_amplifier_bias` (rounded to the nearest
                integer value).
            percentile_cut(float): Percentage of outliers to ignore on both
                sides of the distribution for estimating bias (when bias_subtracted
                is True).
            bias_point(int): When bias_subtracted is True, raw pixel values will
                be offset so that the bias value is shifted to the specified
                value.  Must be between 0 - 0xffff and smaller than the estimated
                bias value.

        Returns:
            numpy.ndarray: 2D array of pixel values as unsigned 16-bit integers.
                The return value is a copy (rather than a view) of the raw data
                array in the file.  When bias_subtracted is True, saturated values
                (0xffff) will remain saturated, and values that underflow (<0) after
                bias subtraction will be set to zero.
        """
        bias_point = int(round(bias_point))
        if bias_point < 0 or bias_point > 0xffff:
            raise ValueError('Invalid 16-bit unsigned bias value: {}.'.format(bias_point))
        regions = []
        for amplifier in range(4):
            data = self.data[self.get_amplifier_region(amplifier,'data')].copy()
            if bias_subtracted:
                overflow = (data == np.uint16(0xffff))
                overscan = self.data[self.get_amplifier_region(amplifier,'overscan')]
                bias_value = self.get_amplifier_bias(amplifier, percentile_cut)
                bias_value = int(round(np.median(overscan)))
                if bias_value < bias_point:
                    raise ValueError('Invalid bias_point ({}) < estimated bias ({}).'
                        .format(bias_point, bias_value))
                offset = np.uint16(bias_value - bias_point)
                underflow = (data < offset)
                data[~overflow & ~underflow] -= offset
                data[underflow] = 0
            regions.append(data)
        return np.hstack(
            (np.vstack((regions[0], regions[2])), np.vstack((regions[1], regions[3]))))

    def read_plug_map(self, speclog_path=None):
        """Read the plug map associated with this plate and exposure.

        Plug maps are not stored under the same SAS tree as other data products
        so reading a plug map requires that a separate directory containing plug
        map files is already available on the local disk. To download the set of
        plug map files used by the pipeline, use the SDSS public SVN repository to
        make a local checkout of the ``speclog`` product::

            svn co https://svn.sdss.org/public/data/sdss/speclog/trunk speclog

        This will take a while to run (about 15 minutes) and will copy about 25 Gb
        of data into the newly created ``speclog`` sub-directory.  Pass the full name
        of this directory to this method.

        You can create your ``speclog`` directory anywhere, but you should avoid
        putting it under \$BOSS_SAS_PATH in your \$BOSS_LOCAL_ROOT since plug map files
        do not originate from SAS. Note that ~11Gb of the files in ``speclog`` are
        `guidermon files
        <http://data.sdss3.org/datamodel/files/SPECLOG_DIR/MJD/guidermon.html>`__
        and not very useful.  You cannot avoid downloading them with the svn
        checkout, but you can delete them after the checkout using::

            cd speclog
            find . -name 'guiderMon-*.par' -delete

        The plug maps are `yanny parameter files
        <https://www.sdss3.org/dr8/software/par.php>`__.
        See the `plug map datamodel
        <http://data.sdss3.org/datamodel/files/PLATELIST_DIR/runs/PLATERUN/plPlugMap.html>`__
        for details.

        Args:
            speclog_path(str): The local path to a directory containing
                plPlugMapM files organized in subdirectories by observation
                MJD.  If None is specified, the value of the BOSS_SPECLOG environment
                variable will be used, if available.

        Returns:
            pydl.pydlutils.yanny.yanny: Object wrapper around the plug map yanny file.
                The PLUGMAPOBJ attribute contains the table of plugging info.
        """
        if speclog_path is None:
            speclog_path = os.getenv('BOSS_SPECLOG')
        if speclog_path is None:
            raise RuntimeError('No speclog_path specified and BOSS_SPECLOG is not set.')

        obs_mjd = '{:d}'.format(self.header['MJD'])
        plug_map_name = 'plPlugMapM-{}.par'.format(self.header['NAME'])
        plug_map_path = os.path.join(speclog_path, obs_mjd, plug_map_name)
        return pydl.pydlutils.yanny.yanny(plug_map_path, np=True)

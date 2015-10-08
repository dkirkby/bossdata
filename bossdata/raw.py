# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Work with raw spectroscopic science and calibration images.
"""

import pydl.pydlutils.yanny


class RawImageFile(object):
    """Wrapper for a raw science or calibration image.

    Each camera (b1,b2,r1,r2) has its own raw image file for a given exposure
    ID. See the `raw image datamodel
    <http://data.sdss3.org/datamodel/files/BOSS_SPECTRO_DATA/MJD/sdR.html>`__
    for details.  The FITS file is opened, read, and closed by the constructor.

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

    def read_plug_map(self, speclog_path):
        """Read the plug map associated with this plate and exposure.

        Plug maps are not stored under the same SAS tree as other data products
        so reading a plug map requires that a separate directory containing plug
        map files is already available on the local disk. To download the set of
        plug map files used by the pipeline, use the SDSS public SVN repository to
        make a local checkout of the ``speclog`` product::

            svn co ​https://svn.sdss.org/public/data/sdss/speclog/trunk speclog

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
                MJD.

        Returns:
            numpy.ndarray: Structured array containing the contents of the PLUGMAPOBJ
                table read from the yanny parameter file.
        """
        obs_mjd = '{:d}'.format(self.header['MJD'])
        plug_map_name = 'plPlugMapM-{}.par'.format(self.header['NAME'])
        plug_map_path = os.path.join(speclog_path, obs_mjd, plug_map_name)
        plug_map = pydl.pydlutils.yanny.yanny(plug_map_path, np=True)
        return plug_map['PLUGMAPOBJ']

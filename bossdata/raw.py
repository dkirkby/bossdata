# -*- coding: utf-8 -*-
# Licensed under a MIT style license - see LICENSE.rst

"""Work with raw spectroscopic science and calibration images.
"""


class RawImageFile(object):
    """Wrapper for a raw science or calibration image.

    Each camera (b1,b2,r1,r2) has its own raw image file for a given exposure
    ID. See the `datamodel
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

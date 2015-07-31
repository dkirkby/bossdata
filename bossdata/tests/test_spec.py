# Licensed under a MIT style license - see LICENSE.rst

import numpy as np

from .. import spec


def test_pixel_index():
    assert spec.get_fiducial_pixel_index(3500.26) == 0
    assert spec.get_fiducial_pixel_index(3564.51) == 79


def test_fiducial_range():
    assert np.array_equal(
        spec.get_fiducial_pixel_index(spec.fiducial_wavelengths),
        np.arange(*spec.fiducial_pixel_index_range))

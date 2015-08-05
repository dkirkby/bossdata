# Licensed under a MIT style license - see LICENSE.rst

import numpy as np

from .. import spec


def test_fiducial_range():
    fiducial_grid = np.power(10., spec.fiducial_loglam)
    assert np.allclose(
        spec.get_fiducial_pixel_index(fiducial_grid),
        np.arange(*spec.fiducial_pixel_index_range))

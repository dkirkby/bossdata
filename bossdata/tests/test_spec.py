# Licensed under a MIT style license - see LICENSE.rst

import pytest

from .. import spec

def test_pixel_index():
    assert spec.get_fiducial_pixel_index(3500.26) == 0
    assert spec.get_fiducial_pixel_index(3564.51) == 79

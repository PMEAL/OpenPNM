import pytest
import numpy as np
import openpnm as op
from numpy.testing import assert_allclose
from openpnm import topotools


class PointsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def teardown_class(self):
        self.ws.clear()

    def test_iscoplanar(self):
        # Generate planar points with several parallel vectors at start
        coords = [[0, 0, 0], [0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, 2]]
        assert topotools.iscoplanar(coords)
        # NON-planar points, also with parallel vectors
        coords = [[0, 0, 0], [0, 0, 0], [0, 0, 1], [0, 0, 2], [1, 1, 2]]
        assert ~topotools.iscoplanar(coords)
        # Planar points, none parallel
        coords = [[0, 0, 0], [0, 1, 2], [0, 2, 1], [0, 3, 2], [0, 2, 3]]
        assert topotools.iscoplanar(coords)
        # Non-planar points, none parallel
        coords = [[0, 0, 0], [0, 1, 2], [0, 2, 1], [0, 3, 3], [1, 1, 2]]
        assert ~topotools.iscoplanar(coords)

    def test_generate_base_points(self):
        pass





if __name__ == '__main__':

    t = PointsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: ' + item)
            t.__getattribute__(item)()

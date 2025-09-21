import numpy as np
from numpy.testing import assert_array_almost_equal

import openpnm as op
import openpnm.models.geometry as gm


class ThroatEndpointsTest:

    def setup_class(self):
        self.net = op.network.Demo(shape=[5, 1, 1], spacing=1)

    def test_spheres_and_cylinders(self):
        self.net.add_model(
            propname='throat.endpoints',
            model=op.models.geometry.throat_endpoints.spheres_and_cylinders
        )
        delta = (self.net['throat.endpoints']['tail']
                 - self.net['throat.endpoints']['head'])[:, 0]
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.circles_and_rectangles)
        assert_array_almost_equal(delta, self.net['throat.length'], decimal=10)


if __name__ == '__main__':

    t = ThroatEndpointsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()

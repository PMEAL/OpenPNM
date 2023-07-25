import numpy as np
from numpy.testing import assert_approx_equal

import openpnm as op
import openpnm.models.geometry.throat_perimeter as mods


class ThroatPerimeterTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.net['throat.diameter'] = 0.1

    def test_cylinder(self):
        self.net.add_model(propname='throat.perimeter',
                           model=mods.cylinder,
                           regen_mode='normal')
        a = np.array([0.31415927])
        b = np.unique(self.net['throat.perimeter'])
        assert_approx_equal(a, b)

    def test_cuboid(self):
        self.net.add_model(propname='throat.perimeter',
                           model=mods.cuboid,
                           regen_mode='normal')
        a = np.array([0.4])
        b = np.unique(self.net['throat.perimeter'])
        assert_approx_equal(a, b)

    def test_rectangle(self):
        self.net.add_model(propname='throat.perimeter',
                           model=mods.rectangle,
                           regen_mode='normal')
        a = np.array([1.0])
        b = np.unique(self.net['throat.perimeter'])
        assert_approx_equal(a, b)


if __name__ == '__main__':

    t = ThroatPerimeterTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()

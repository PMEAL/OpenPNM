import numpy as np
from numpy.testing import assert_approx_equal

import openpnm as op
import openpnm.models.geometry.throat_cross_sectional_area as mods


class ThroatCrossSectionalAreaTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.net['pore.diameter'] = 0.2
        self.net['throat.length'] = 1
        self.net['throat.diameter'] = 0.1

    def test_sphere(self):
        self.net.add_model(propname='throat.cross_sectional_area',
                           model=mods.cylinder,
                           regen_mode='normal')
        a = np.array([0.007853981])
        b = np.unique(self.net['throat.cross_sectional_area'])
        assert_approx_equal(a, b)

    def test_cube(self):
        self.net.add_model(propname='throat.cross_sectional_area',
                              model=mods.cuboid,
                              regen_mode='normal')
        a = np.array([0.01])
        b = np.unique(self.net['throat.cross_sectional_area'])
        assert_approx_equal(a, b)

    def test_rectangle(self):
        self.net.add_model(propname='throat.cross_sectional_area',
                           model=mods.rectangle,
                           regen_mode='normal')
        a = np.array([0.1])
        b = np.unique(self.net['throat.cross_sectional_area'])
        assert_approx_equal(a, b)


if __name__ == '__main__':

    t = ThroatCrossSectionalAreaTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()

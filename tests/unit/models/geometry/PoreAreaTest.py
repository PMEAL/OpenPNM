import openpnm as op
import openpnm.models.geometry.pore_cross_sectional_area as mods
import numpy as np
from numpy.testing import assert_approx_equal


class PoreAreaTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        self.geo['pore.diameter'] = 1.0
        self.geo['throat.area'] = 0.1

    def test_sphere(self):
        self.geo.add_model(propname='pore.area',
                           model=mods.sphere,
                           regen_mode='normal')
        a = np.array([0.78539816])
        b = np.unique(self.geo['pore.area'])
        assert_approx_equal(a, b)

    def test_cube(self):
        self.geo.add_model(propname='pore.area',
                           model=mods.cube,
                           regen_mode='normal')
        a = np.array([1.0])
        b = np.unique(self.geo['pore.area'])
        assert_approx_equal(a, b)

    def test_circle(self):
        self.geo.add_model(propname='pore.area',
                           model=mods.circle,
                           regen_mode='normal')
        a = np.array([1.0])
        b = np.unique(self.geo['pore.area'])
        assert_approx_equal(a, b)

    def test_square(self):
        self.geo.add_model(propname='pore.area',
                           model=mods.square,
                           regen_mode='normal')
        a = np.array([1.0])
        b = np.unique(self.geo['pore.area'])
        assert_approx_equal(a, b)


if __name__ == '__main__':

    t = PoreAreaTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

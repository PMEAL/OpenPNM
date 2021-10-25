import openpnm as op
import openpnm.models.geometry.throat_cross_sectional_area as mods
import numpy as np
from numpy.testing import assert_approx_equal


class ThroatSurfaceAreaTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        self.geo['pore.diameter'] = 0.2
        self.geo['throat.length'] = 1
        self.geo['throat.diameter'] = 0.1

    def test_sphere(self):
        self.geo.add_model(propname='throat.area',
                           model=mods.cylinder,
                           regen_mode='normal')
        a = np.array([0.007853981])
        b = np.unique(self.geo['throat.area'])
        assert_approx_equal(a, b)

    def test_cube(self):
        self.geo.add_model(propname='throat.area',
                           model=mods.cuboid,
                           regen_mode='normal')
        a = np.array([0.01])
        b = np.unique(self.geo['throat.area'])
        assert_approx_equal(a, b)

    def test_rectangle(self):
        self.geo.add_model(propname='throat.area',
                           model=mods.rectangle,
                           regen_mode='normal')
        a = np.array([0.1])
        b = np.unique(self.geo['throat.area'])
        assert_approx_equal(a, b)


if __name__ == '__main__':

    t = ThroatSurfaceAreaTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

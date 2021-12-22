import openpnm as op
import openpnm.models.geometry.throat_perimeter as mods
import numpy as np
from numpy.testing import assert_approx_equal


class ThroatPerimeterTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.air = op.phase.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        self.geo['throat.diameter'] = 0.1

    def test_cylinder(self):
        self.geo.add_model(propname='throat.perimeter',
                           model=mods.cylinder,
                           regen_mode='normal')
        a = np.array([0.31415927])
        b = np.unique(self.geo['throat.perimeter'])
        assert_approx_equal(a, b)

    def test_cuboid(self):
        self.geo.add_model(propname='throat.perimeter',
                           model=mods.cuboid,
                           regen_mode='normal')
        a = np.array([0.4])
        b = np.unique(self.geo['throat.perimeter'])
        assert_approx_equal(a, b)

    def test_rectangle(self):
        self.geo.add_model(propname='throat.perimeter',
                           model=mods.rectangle,
                           regen_mode='normal')
        a = np.array([1.0])
        b = np.unique(self.geo['throat.perimeter'])
        assert_approx_equal(a, b)


if __name__ == '__main__':

    t = ThroatPerimeterTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

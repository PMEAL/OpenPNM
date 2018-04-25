import openpnm as op
import openpnm.models.geometry.throat_length as mods
import numpy as np
from numpy.testing import assert_approx_equal


class ThroatStraightTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        self.geo['pore.diameter'] = .05

    def test_straight(self):
        self.geo.add_model(propname='throat.length',
                           model=mods.straight,
                           regen_mode='normal')
        a = np.array([0.95])
        b = np.unique(self.geo['throat.length'])
        assert_approx_equal(a, b)


if __name__ == '__main__':

    t = ThroatStraightTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

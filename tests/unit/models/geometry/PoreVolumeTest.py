import openpnm as op
import openpnm.models.geometry.pore_volume as mods
import numpy as np
from numpy.testing import assert_approx_equal


class PoreVolumeTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.air = op.phase.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        self.geo['pore.diameter'] = 1.05
        self.geo['throat.area'] = 0.1

    def test_sphere(self):
        self.geo.add_model(propname='pore.volume',
                           model=mods.sphere,
                           regen_mode='normal')
        a = np.array([0.52359878*1.05**3])
        b = np.unique(self.geo['pore.volume'])
        assert_approx_equal(a, b)

    def test_circle(self):
        self.geo.add_model(propname='pore.volume',
                           model=mods.circle,
                           regen_mode='normal')
        a = np.array([3.14159265/4*1.05**2])
        b = np.unique(self.geo['pore.volume'])
        assert_approx_equal(a, b)

    def test_cube(self):
        self.geo.add_model(propname='pore.volume',
                           model=mods.cube,
                           regen_mode='normal')
        a = np.array([1.0*1.05**3])
        b = np.unique(self.geo['pore.volume'])
        assert_approx_equal(a, b)

    def test_square(self):
        self.geo.add_model(propname='pore.volume',
                           model=mods.square,
                           regen_mode='normal')
        a = np.array([1.0*1.05**2])
        b = np.unique(self.geo['pore.volume'])
        assert_approx_equal(a, b)

    def test_effective(self):
        net = op.network.Cubic(shape=[2, 2, 1])
        net['pore.volume'] = 0.5
        net['throat.volume'] = 0.25
        net.add_model(propname='pore.volume_effective',
                      model=mods.effective)
        a = np.array([0.5 + 0.25])
        b = np.unique(net['pore.volume_effective'])
        assert_approx_equal(a, b)


if __name__ == '__main__':

    t = PoreVolumeTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

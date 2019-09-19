import openpnm as op
import openpnm.models.geometry.throat_volume as mods
import numpy as np
from numpy.testing import assert_approx_equal


class ThroatVolumeTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        self.geo['throat.diameter'] = 0.1
        self.geo['throat.length'] = 1.0
        self.geo['throat.area'] = 0.03

    def test_lens_and_pendular_ring(self):
        net = op.network.Cubic(shape=[2, 1, 1])
        net['pore.diameter'] = 0.5
        net['throat.diameter'] = 0.25
        mod = op.models.geometry.throat_volume.lens
        net.add_model(propname='throat.lens_volume',
                      model=mod)
        Vlens = net['throat.lens_volume']
        assert Vlens == 2*(0.006733852203712552)
        mod = op.models.geometry.throat_volume.pendular_ring
        net.add_model(propname='throat.ring_volume',
                      model=mod)
        Vcyl = 2*(0.01315292522620208)
        Vring = net['throat.ring_volume']
        assert (Vcyl - Vring) == 2*(0.006733852203712552)

    def test_cylinder(self):
        self.geo.add_model(propname='throat.volume',
                           model=mods.cylinder,
                           regen_mode='normal')
        a = np.array([0.007853981])
        b = np.unique(self.geo['throat.volume'])
        assert_approx_equal(a, b)

    def test_cube(self):
        self.geo.add_model(propname='throat.volume',
                           model=mods.cuboid)
        a = np.array([0.01])
        b = np.unique(self.geo['throat.volume'])
        assert_approx_equal(a, b)

    def test_rectangle(self):
        self.geo.add_model(propname='throat.volume',
                           model=mods.rectangle)
        a = np.array([0.1])
        b = np.unique(self.geo['throat.volume'])
        assert_approx_equal(a, b)

    def test_extrusion(self):
        self.geo.add_model(propname='throat.volume',
                           throat_area='throat.area',
                           model=mods.extrusion)
        a = np.array([0.03])
        b = np.unique(self.geo['throat.volume'])
        assert np.allclose(a, b)


if __name__ == '__main__':

    t = ThroatVolumeTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

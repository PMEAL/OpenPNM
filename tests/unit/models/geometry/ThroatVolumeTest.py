import openpnm as op
import openpnm.models.geometry.throat_volume as mods
import numpy as np
from numpy.testing import assert_approx_equal


class ThroatVolumeTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.net['throat.diameter'] = 0.1
        self.net['throat.length'] = 1.0
        self.net['throat.cross_sectional_area'] = 0.03

    def test_lens_and_pendular_ring(self):
        net = op.network.Cubic(shape=[2, 1, 1])
        net['pore.diameter'] = 0.5
        net['throat.diameter'] = 0.25
        mod = op.models.geometry.throat_volume.lens
        net.add_model(propname='throat.lens_volume',
                      model=mod)
        Vlens = net['throat.lens_volume']
        assert np.isclose(Vlens, 2*0.00084173)
        mod = op.models.geometry.throat_volume.pendular_ring
        net.add_model(propname='throat.ring_volume',
                      model=mod)
        Vcyl = 2*(0.00164412)
        Vring = net['throat.ring_volume']
        assert np.isclose(Vcyl - Vring, 2*0.00084173)

    def test_cylinder(self):
        self.net.add_model(propname='throat.volume',
                           model=mods.cylinder,
                           regen_mode='normal')
        a = np.array([0.007853981])
        b = np.unique(self.net['throat.volume'])
        assert_approx_equal(a, b)

    def test_cube(self):
        self.net.add_model(propname='throat.volume',
                           model=mods.cuboid)
        a = np.array([0.01])
        b = np.unique(self.net['throat.volume'])
        assert_approx_equal(a, b)

    def test_rectangle(self):
        self.net.add_model(propname='throat.volume',
                           model=mods.rectangle)
        a = np.array([0.1])
        b = np.unique(self.net['throat.volume'])
        assert_approx_equal(a, b)

    def test_extrusion(self):
        self.net.add_model(propname='throat.volume',
                           throat_area='throat.cross_sectional_area',
                           model=mods.extrusion)
        a = np.array([0.03])
        b = np.unique(self.net['throat.volume'])
        assert np.allclose(a, b)


if __name__ == '__main__':

    t = ThroatVolumeTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

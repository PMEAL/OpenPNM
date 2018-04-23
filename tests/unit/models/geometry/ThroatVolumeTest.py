import openpnm as op
import openpnm.models.geometry.throat_volume as mods
import numpy as np


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

    def test_cylinder(self):
        self.geo.add_model(propname='throat.volume',
                           model=mods.cylinder,
                           regen_mode='eager')
        a = np.array([0.00785398])
        b = np.unique(self.geo['throat.volume'])
        assert np.allclose(a, b)

    def test_cube(self):
        self.geo.add_model(propname='throat.volume',
                           model=mods.cuboid,
                           regen_mode='eager')
        a = np.array([0.01])
        b = np.unique(self.geo['throat.volume'])
        assert np.allclose(a, b)

    def test_extrusion(self):
        self.geo.add_model(propname='throat.volume',
                           throat_area='throat.area',
                           model=mods.extrusion,
                           regen_mode='eager')
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

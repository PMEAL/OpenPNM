import openpnm as op
import openpnm.models.geometry.pore_surface_area as mods
import numpy as np


class PoreSurfaceAreaTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        self.geo['pore.diameter'] = 1
        self.geo['throat.area'] = 0.1

    def test_sphere(self):
        self.geo.add_model(propname='pore.area',
                           model=mods.sphere,
                           regen_mode='normal')
        a = np.array([2.54159265, 2.64159265, 2.74159265, 2.84159265])
        b = np.unique(self.geo['pore.area'])
        assert np.allclose(a, b)

    def test_cube(self):
        self.geo.add_model(propname='pore.area',
                           model=mods.cube,
                           regen_mode='normal')
        a = np.array([5.4, 5.5, 5.6, 5.7])
        b = np.unique(self.geo['pore.area'])
        assert np.allclose(a, b)


if __name__ == '__main__':

    t = PoreSurfaceAreaTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

import openpnm as op
import openpnm.models.geometry.pore_surface_area as mods
import numpy as np
from numpy.testing import assert_allclose


class PoreSurfaceAreaTest:
    def setup_class(self):
        self.net3D = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.net3D['pore.diameter'] = 1
        self.net3D['throat.cross_sectional_area'] = 0.1
        self.net2D = op.network.Cubic(shape=[5, 5, 1], spacing=1.0)
        self.net2D['pore.diameter'] = 1
        self.net2D['throat.cross_sectional_area'] = 0.1

    def test_sphere(self):
        self.net3D.add_model(propname='pore.surface_area',
                             model=mods.sphere,
                             regen_mode='normal')
        a = np.array([2.54159265, 2.64159265, 2.74159265, 2.84159265])
        b = np.unique(self.net3D['pore.surface_area'])
        assert_allclose(a, b)

    def test_circle(self):
        self.net2D.add_model(propname='pore.surface_area',
                             model=mods.circle,
                             regen_mode='normal')
        a = np.array([2.74159265, 2.84159265, 2.94159265])
        b = np.unique(self.net2D['pore.surface_area'])
        assert_allclose(a, b)

    def test_cube(self):
        self.net3D.add_model(propname='pore.surface_area',
                             model=mods.cube,
                             regen_mode='normal')
        a = np.array([5.4, 5.5, 5.6, 5.7])
        b = np.unique(self.net3D['pore.surface_area'])
        assert_allclose(a, b)

    def test_square(self):
        self.net2D.add_model(propname='pore.surface_area',
                             model=mods.square,
                             regen_mode='normal')
        a = np.array([3.6, 3.7, 3.8])
        b = np.unique(self.net2D['pore.surface_area'])
        assert_allclose(a, b)

    def test_circle_multi_geom(self):
        net = op.network.Cubic(shape=[10, 1, 1], spacing=1.0)
        net['pore.left'] = False
        net['pore.left'][[0, 1, 2]] = True
        net['pore.diameter'] = 1
        net['throat.cross_sectional_area'] = 0.1
        net['throat.cross_sectional_area'][[0, 1, 2]] = 0.3
        net['pore.right'] = ~net['pore.left']
        net.add_model(propname='pore.surface_area',
                      model=mods.circle,
                      domain='left',
                      regen_mode='normal')
        net.add_model(propname='pore.surface_area',
                      model=mods.circle,
                      domain='right',
                      regen_mode='normal')
        a = np.array([2.84159265, 2.54159265, 2.54159265])
        b = net['pore.surface_area@left']
        c = np.array([2.74159265, 2.94159265, 2.94159265, 2.94159265,
                      2.94159265, 2.94159265, 3.04159265])
        d = net['pore.surface_area@right']
        assert_allclose(a, b)
        assert_allclose(c, d)


if __name__ == '__main__':

    t = PoreSurfaceAreaTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()

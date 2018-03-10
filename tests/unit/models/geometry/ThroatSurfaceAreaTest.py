import openpnm as op
import scipy as sp


class ThroatSurfaceAreaTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['throat.diameter'] = sp.ones((self.geo.Nt, ))
        self.geo['throat.length'] = sp.ones((self.geo.Nt, ))
        self.geo['throat.perimeter'] = sp.ones((self.geo.Nt, ))

    def test_cylinder(self):
        f = op.models.geometry.throat_surface_area.cylinder
        self.geo.add_model(propname='throat.surface_area',
                           model=f)
        self.geo.regenerate_models()
        assert sp.all(self.geo['throat.surface_area'] == sp.pi)

    def test_cuboid(self):
        f = op.models.geometry.throat_surface_area.cuboid
        self.geo.add_model(propname='throat.surface_area',
                           model=f,)
        self.geo.regenerate_models()
        assert sp.all(self.geo['throat.surface_area'] == 4)

    def test_extrusion(self):
        f = op.models.geometry.throat_surface_area.extrusion
        self.geo.add_model(propname='throat.surface_area',
                           model=f)
        self.geo.regenerate_models()
        assert sp.all(self.geo['throat.surface_area'] == 1)


if __name__ == '__main__':

    t = ThroatSurfaceAreaTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

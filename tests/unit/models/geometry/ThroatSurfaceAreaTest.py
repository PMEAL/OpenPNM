import OpenPNM
import scipy as sp


class ThroatSurfaceAreaTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts)
        self.geo['throat.diameter'] = sp.ones((self.geo.Nt, ))
        self.geo['throat.length'] = sp.ones((self.geo.Nt, ))
        self.geo['throat.perimeter'] = sp.ones((self.geo.Nt, ))

    def test_cylinder(self):
        f = OpenPNM.Geometry.models.throat_surface_area.cylinder
        self.geo.models.add(propname='throat.surface_area',
                            model=f)
        assert sp.all(self.geo['throat.surface_area'] == sp.pi)

    def test_cuboid(self):
        f = OpenPNM.Geometry.models.throat_surface_area.cuboid
        self.geo.models.add(propname='throat.surface_area',
                            model=f,)
        assert sp.all(self.geo['throat.surface_area'] == 4)

    def test_extrusion(self):
        f = OpenPNM.Geometry.models.throat_surface_area.extrusion
        self.geo.models.add(propname='throat.surface_area',
                            model=f)
        assert sp.all(self.geo['throat.surface_area'] == 1)


if __name__ == '__main__':

    t = ThroatSizeTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
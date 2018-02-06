import OpenPNM
import scipy as sp
import scipy.stats as spst
import OpenPNM.Geometry.models.pore_diameter as mods


class PoreDiameterTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts)
        self.geo['pore.seed'] = sp.rand(self.geo.Np)

    def test_normal(self):
        self.geo.models.add(propname='pore.diameter',
                            model=mods.normal,
                            scale=0.01,
                            loc=0.5,
                            seeds='pore.seed')
        assert 0.45 < sp.mean(self.geo['pore.diameter']) < 0.55
        del self.geo['pore.diameter']

    def test_weibull(self):
        self.geo.models.add(propname='pore.diameter',
                            model=mods.weibull,
                            shape=1.5,
                            scale=0.0001,
                            loc=0.001,
                            seeds='pore.seed')
        assert sp.amin(self.geo['pore.diameter']) > 0.001
        del self.geo['pore.diameter']

    def test_generic(self):
        func = spst.gamma(a=2, loc=0.001, scale=0.0001)
        self.geo.models.add(propname='pore.diameter',
                            model=mods.generic,
                            func=func,
                            seeds='pore.seed')
        assert sp.amin(self.geo['pore.diameter']) > 0.001
        del self.geo['pore.diameter']

    def test_largest_sphere(self):
        net = OpenPNM.Network.Cubic(shape=[5, 5, 5], spacing=[0.1, 0.2, 0.3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net, pores=net.Ps,
                                               throats=net.Ts)
        geo.models.add(propname='pore.diameter',
                       model=mods.largest_sphere,
                       iters=1)
        dmin = sp.amin(geo['pore.diameter'])
        assert dmin <= 0.1
        geo.models['pore.diameter']['iters'] = 5
        geo.regenerate()
        assert dmin < sp.amin(geo['pore.diameter'])

    def test_largest_sphere_multiple_geometries(self):
        net = OpenPNM.Network.Cubic(shape=[10, 10, 10], spacing=[5, 5, 5])
        net['pore.coords'][net.pores('top')] += [0, 0, -3]
        geom2 = OpenPNM.Geometry.GenericGeometry(network=net,
                                                 pores=net.pores('top'))
        geom2['pore.diameter'] = 1.0
        Ps = net.pores('top', mode='not')
        geom1 = OpenPNM.Geometry.GenericGeometry(network=net,
                                                 pores=Ps,
                                                 throats=net.Ts)
        mod = OpenPNM.Geometry.models.pore_diameter.largest_sphere
        geom1.models.add(propname='pore.diameter',
                         model=mod,
                         iters=15)
        assert sp.all(geom2['pore.diameter'] == 1.0)
        assert sp.all(sp.ceil(sp.unique(geom1['pore.diameter'])) == [3.0, 5.0])

    def test_equivalent_cube(self):
        net = OpenPNM.Network.Cubic(shape=[5, 5, 5], spacing=[0.1, 0.2, 0.3])
        geo = OpenPNM.Geometry.GenericGeometry(network=net, pores=net.Ps,
                                               throats=net.Ts)
        geo['pore.volume'] = 8.0
        geo.models.add(propname='pore.diameter',
                       model=mods.equivalent_cube)
        assert sp.all(geo['pore.diameter'] == 2.0)

if __name__ == '__main__':
    t = PoreDiameterTest()
    t.setup_class()
    t.test_equivalent_cube()

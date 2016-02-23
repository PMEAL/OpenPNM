import OpenPNM
import scipy as sp
import scipy.stats as spst


class PoreDiameterTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts)
        self.geo['pore.seed'] = sp.rand(self.geo.Np)

    def test_normal(self):
        import OpenPNM.Geometry.models.pore_diameter as mods
        self.geo.models.add(propname='pore.diameter',
                            model=mods.normal,
                            scale=0.01,
                            loc=0.5,
                            seeds='pore.seed')
        assert 0.45 < sp.mean(self.geo['pore.diameter']) < 0.55
        del self.geo['pore.diameter']

    def test_weibull(self):
        import OpenPNM.Geometry.models.pore_diameter as mods
        self.geo.models.add(propname='pore.diameter',
                            model=mods.weibull,
                            shape=1.5,
                            scale=0.0001,
                            loc=0.001,
                            seeds='pore.seed')
        assert sp.amin(self.geo['pore.diameter']) > 0.001
        del self.geo['pore.diameter']

    def test_generic(self):
        import OpenPNM.Geometry.models.pore_diameter as mods
        func = spst.gamma(a=2, loc=0.001, scale=0.0001)
        self.geo.models.add(propname='pore.diameter',
                            model=mods.generic,
                            func=func,
                            seeds='pore.seed')
        assert sp.amin(self.geo['pore.diameter']) > 0.001
        del self.geo['pore.diameter']

import OpenPNM
from OpenPNM.Geometry import models as gm
import scipy as sp
import scipy.stats as spst


class ThroatDiameterTest:

    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[2, 2, 2])
        self.net.add_boundaries()
        Ps = self.net.pores('boundary', mode='not')
        Ts = self.net.throats('boundary', mode='not')
        self.geo = OpenPNM.Geometry.TestGeometry(network=self.net,
                                                 pores=Ps,
                                                 throats=Ts)
        BPs = self.net.pores('boundary')
        BTs = self.net.throats('boundary')
        self.boun = OpenPNM.Geometry.Boundary(network=self.net,
                                              pores=BPs,
                                              throats=BTs)

    def test_min_pore(self):
        self.geo.models.add(propname='throat.diameter',
                            model=gm.throat_diameter.minpore)
        self.boun.models.add(propname='throat.diameter',
                             model=gm.throat_diameter.minpore)
        BPs = self.net["throat.conns"][self.net.throats('boundary')][:, 0]
        assert sp.sum(self.net["pore.diameter"][BPs] -
                      self.boun["throat.diameter"]) == 0.0

    def test_normal(self):
        import OpenPNM.Geometry.models.pore_diameter as mods
        self.geo.models.add(propname='throat.diameter',
                            model=mods.normal,
                            scale=0.01,
                            loc=0.5,
                            seeds='throat.seed')
        assert 0.45 < sp.mean(self.geo['throat.diameter']) < 0.55
        del self.geo['throat.diameter']

    def test_weibull(self):
        import OpenPNM.Geometry.models.pore_diameter as mods
        self.geo.models.add(propname='throat.diameter',
                            model=mods.weibull,
                            shape=1.5,
                            scale=0.0001,
                            loc=0.001,
                            seeds='throat.seed')
        assert sp.amin(self.geo['throat.diameter']) > 0.001
        del self.geo['throat.diameter']

    def test_generic(self):
        import OpenPNM.Geometry.models.pore_diameter as mods
        func = spst.gamma(a=2, loc=0.001, scale=0.0001)
        self.geo.models.add(propname='throat.diameter',
                            model=mods.generic,
                            func=func,
                            seeds='throat.seed')
        assert sp.amin(self.geo['throat.diameter']) > 0.001
        del self.geo['throat.diameter']

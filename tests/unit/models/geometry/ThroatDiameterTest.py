import openpnm as op
import openpnm.models.geometry as gm
import scipy as sp
import scipy.stats as spst


class ThroatSizeTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[2, 2, 2])
        self.net.add_boundary_pores()
        Ps = self.net.pores('*boundary', mode='not')
        Ts = self.net.throats('*boundary', mode='not')
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=Ps,
                                               throats=Ts)
        sp.random.RandomState(seed=0)
        self.geo['throat.seed'] = sp.rand(self.geo.Nt)
        BPs = self.net.pores('*boundary')
        BTs = self.net.throats('*boundary')
        self.boun = op.geometry.GenericGeometry(network=self.net,
                                                pores=BPs,
                                                throats=BTs)

    def test_normal(self):
        self.geo.add_model(propname='throat.diameter',
                           model=gm.throat_size.normal,
                           scale=0.01,
                           loc=0.5,
                           seeds='throat.seed')
        assert 0.45 < sp.mean(self.geo['throat.diameter']) < 0.55
        del self.geo['throat.diameter']

    def test_weibull(self):
        self.geo.add_model(propname='throat.diameter',
                           model=gm.throat_size.weibull,
                           shape=1.5,
                           scale=0.0001,
                           loc=0.001,
                           seeds='throat.seed')
        assert sp.amin(self.geo['throat.diameter']) > 0.001
        del self.geo['throat.diameter']

    def test_generic(self):
        func = spst.gamma(a=2, loc=0.001, scale=0.0001)
        self.geo.add_model(propname='throat.diameter',
                           model=gm.throat_size.generic,
                           func=func,
                           seeds='throat.seed')
        assert sp.amin(self.geo['throat.diameter']) > 0.001
        del self.geo['throat.diameter']


if __name__ == '__main__':

    t = ThroatSizeTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

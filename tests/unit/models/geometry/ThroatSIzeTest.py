import numpy as np
import openpnm as op
import scipy.stats as spst
import openpnm.models.geometry as gm


class ThroatSizeTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[2, 2, 2])
        self.net.add_boundary_pores()
        Ps = self.net.pores('*boundary', mode='nor')
        Ts = self.net.throats('*boundary', mode='nor')
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=Ps,
                                               throats=Ts)
        np.random.RandomState(seed=0)
        self.geo['throat.seed'] = np.random.rand(self.geo.Nt)
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
        assert 0.45 < np.mean(self.geo['throat.diameter']) < 0.55
        del self.geo['throat.diameter']

    def test_weibull(self):
        self.geo.add_model(propname='throat.diameter',
                           model=gm.throat_size.weibull,
                           shape=1.5,
                           scale=0.0001,
                           loc=0.001,
                           seeds='throat.seed')
        assert np.amin(self.geo['throat.diameter']) > 0.001
        del self.geo['throat.diameter']

    def test_generic_distribution(self):
        func = spst.gamma(a=2, loc=0.001, scale=0.0001)
        self.geo.add_model(propname='throat.diameter',
                           model=gm.throat_size.generic_distribution,
                           func=func,
                           seeds='throat.seed')
        assert np.amin(self.geo['throat.diameter']) > 0.001
        del self.geo['throat.diameter']

    def test_from_neighbor_pores(self):
        self.geo['pore.diameter'] = 0.1
        self.geo.add_model(propname='throat.diameter',
                           model=gm.throat_size.from_neighbor_pores,
                           prop='pore.diameter')
        a = np.unique(self.geo['throat.diameter'])
        b = np.array(0.1, ndmin=1)
        assert np.allclose(a, b)
        del self.geo['throat.diameter'], self.geo.models['throat.diameter']

    def test_equivalent_diameter(self):
        self.geo['throat.area'] = 1.0
        self.geo.add_model(propname='throat.diameter',
                           model=gm.throat_size.equivalent_diameter,
                           throat_area='throat.area',
                           throat_shape='circle')
        a = np.unique(self.geo['throat.diameter'])
        b = np.array(1.12837917, ndmin=1)
        assert np.allclose(a, b)
        del self.geo['throat.diameter'], self.geo.models['throat.diameter']
        self.geo.add_model(propname='throat.diameter',
                           model=gm.throat_size.equivalent_diameter,
                           throat_area='throat.area',
                           throat_shape='square')
        a = np.unique(self.geo['throat.diameter'])
        b = np.array(1.0, ndmin=1)
        assert np.allclose(a, b)
        del self.geo['throat.diameter'], self.geo.models['throat.diameter']


if __name__ == '__main__':

    t = ThroatSizeTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

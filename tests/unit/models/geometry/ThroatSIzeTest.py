import numpy as np
import openpnm as op
import scipy.stats as spst
import openpnm.models.geometry as gm


class ThroatSizeTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[2, 2, 2])
        np.random.RandomState(seed=0)
        self.net['throat.seed'] = np.random.rand(self.net.Nt)

    def test_normal(self):
        self.net.add_model(propname='throat.diameter',
                           model=gm.throat_size.normal,
                           scale=0.01,
                           loc=0.5,
                           seeds='throat.seed')
        assert 0.45 < np.mean(self.net['throat.diameter']) < 0.55
        del self.net['throat.diameter']

    def test_weibull(self):
        self.net.add_model(propname='throat.diameter',
                           model=gm.throat_size.weibull,
                           shape=1.5,
                           scale=0.0001,
                           loc=0.001,
                           seeds='throat.seed')
        assert np.amin(self.net['throat.diameter']) > 0.001
        del self.net['throat.diameter']

    def test_generic_distribution(self):
        func = spst.gamma(a=2, loc=0.001, scale=0.0001)
        self.net.add_model(propname='throat.diameter',
                           model=gm.throat_size.generic_distribution,
                           func=func,
                           seeds='throat.seed')
        assert np.amin(self.net['throat.diameter']) > 0.001
        del self.net['throat.diameter']

    def test_from_neighbor_pores(self):
        self.net['pore.diameter'] = 0.1
        self.net.add_model(propname='throat.diameter',
                           model=gm.throat_size.from_neighbor_pores,
                           prop='pore.diameter',
                           domain='all')
        a = np.unique(self.net['throat.diameter'])
        b = np.array(0.1, ndmin=1)
        assert np.allclose(a, b)
        del self.net['throat.diameter'], self.net.models['throat.diameter@all']

    def test_equivalent_diameter(self):
        self.net['throat.cross_sectional_area'] = 1.0
        self.net.add_model(propname='throat.diameter',
                           model=gm.throat_size.equivalent_diameter,
                           throat_area='throat.cross_sectional_area',
                           throat_shape='circle',
                           domain='all')
        a = np.unique(self.net['throat.diameter'])
        b = np.array(1.12837917, ndmin=1)
        assert np.allclose(a, b)
        del self.net['throat.diameter'], self.net.models['throat.diameter@all']
        self.net.add_model(propname='throat.diameter',
                           model=gm.throat_size.equivalent_diameter,
                           throat_area='throat.cross_sectional_area',
                           throat_shape='square',
                           domain='all')
        a = np.unique(self.net['throat.diameter'])
        b = np.array(1.0, ndmin=1)
        assert np.allclose(a, b)
        del self.net['throat.diameter'], self.net.models['throat.diameter@all']


if __name__ == '__main__':

    t = ThroatSizeTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

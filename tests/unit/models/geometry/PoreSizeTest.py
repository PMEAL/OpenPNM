import numpy as np
import openpnm as op
import scipy.stats as spst
import openpnm.models.geometry.pore_size as mods


class PoreSizeTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.net['pore.seed'] = np.random.rand(self.net.Np)

    def test_normal(self):
        self.net.add_model(propname='pore.diameter',
                           model=mods.normal,
                           scale=0.01,
                           loc=0.5,
                           seeds='pore.seed')
        assert 0.45 < np.mean(self.net['pore.diameter']) < 0.55
        del self.net['pore.diameter']

    def test_weibull(self):
        self.net.add_model(propname='pore.diameter',
                           model=mods.weibull,
                           shape=1.5,
                           scale=0.0001,
                           loc=0.001,
                           seeds='pore.seed')
        assert np.amin(self.net['pore.diameter']) > 0.001
        del self.net['pore.diameter']

    def test_generic_distribution(self):
        func = spst.gamma(a=2, loc=0.001, scale=0.0001)
        self.net.add_model(propname='pore.diameter',
                           model=mods.generic_distribution,
                           func=func,
                           seeds='pore.seed')
        assert np.amin(self.net['pore.diameter']) > 0.001
        del self.net['pore.diameter']

    def test_largest_sphere(self):
        net = op.network.Cubic(shape=[5, 5, 5], spacing=[0.1, 0.2, 0.3])
        net.add_model(propname='pore.diameter',
                      model=mods.largest_sphere,
                      iters=5,
                      domain='all')
        dmin = np.amin(net['pore.diameter'])
        assert dmin <= 0.1
        net.models['pore.diameter@all']['iters'] = 5
        net.regenerate_models()
        assert dmin <= np.amin(net['pore.diameter'])

    def test_largest_sphere_multiple_geometries(self):
        net = op.network.Cubic(shape=[5, 5, 5], spacing=[5, 5, 5])
        net['pore.coords'][net.pores('top')] += [0, 0, -3.0]
        net['pore.fixed_diameter@top'] = 1
        net['pore.not_top'] = False
        net['pore.not_top'][net.pores('top', mode='nor')] = True
        net.add_model(propname='pore.diameter',
                      model=mods.largest_sphere,
                      domain='not_top',
                      iters=0)
        assert np.all(np.unique(net['pore.diameter@not_top']) == [1.5, 5.0])
        net['pore.fixed_diameter@top'] = 6
        net.regenerate_models()
        assert np.amin(net['pore.diameter@not_top']) < 0

    def test_equivalent_diameter(self):
        mod = op.models.geometry.pore_size.equivalent_diameter
        self.net['pore.volume'] = 1.0
        self.net.add_model(propname='pore.diameter',
                           model=mod,
                           pore_volume='pore.volume',
                           pore_shape='sphere',
                           domain='all')
        a = np.unique(self.net['pore.diameter'])
        b = np.array(1.24070098, ndmin=1)
        assert np.allclose(a, b)
        del self.net['pore.diameter'], self.net.models['pore.diameter@all']
        self.net.add_model(propname='pore.diameter',
                           model=mod,
                           pore_volume='pore.volume',
                           pore_shape='cube')
        a = np.unique(self.net['pore.diameter'])
        b = np.array(1.0, ndmin=1)
        assert np.allclose(a, b)


if __name__ == '__main__':
    import py
    t = PoreSizeTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())

import py
import numpy as np
import openpnm as op
import scipy.stats as spst
import openpnm.models.geometry.pore_size as mods


class PoreSizeTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.seed'] = np.random.rand(self.geo.Np)

    def test_normal(self):
        self.geo.add_model(propname='pore.diameter',
                           model=mods.normal,
                           scale=0.01,
                           loc=0.5,
                           seeds='pore.seed')
        assert 0.45 < np.mean(self.geo['pore.diameter']) < 0.55
        del self.geo['pore.diameter']

    def test_weibull(self):
        self.geo.add_model(propname='pore.diameter',
                           model=mods.weibull,
                           shape=1.5,
                           scale=0.0001,
                           loc=0.001,
                           seeds='pore.seed')
        assert np.amin(self.geo['pore.diameter']) > 0.001
        del self.geo['pore.diameter']

    def test_generic_distribution(self):
        func = spst.gamma(a=2, loc=0.001, scale=0.0001)
        self.geo.add_model(propname='pore.diameter',
                           model=mods.generic_distribution,
                           func=func,
                           seeds='pore.seed')
        assert np.amin(self.geo['pore.diameter']) > 0.001
        del self.geo['pore.diameter']

    def test_largest_sphere(self):
        net = op.network.Cubic(shape=[5, 5, 5], spacing=[0.1, 0.2, 0.3])
        geo = op.geometry.GenericGeometry(network=net, pores=net.Ps,
                                          throats=net.Ts)
        geo.add_model(propname='pore.diameter',
                      model=mods.largest_sphere,
                      iters=5)
        dmin = np.amin(geo['pore.diameter'])
        assert dmin <= 0.1
        geo.models['pore.diameter']['iters'] = 5
        geo.regenerate_models()
        assert dmin <= np.amin(geo['pore.diameter'])

    def test_largest_sphere_multiple_geometries(self):
        net = op.network.Cubic(shape=[5, 5, 5], spacing=[5, 5, 5])
        net['pore.coords'][net.pores('top')] += [0, 0, -3.0]
        geom2 = op.geometry.GenericGeometry(network=net,
                                            pores=net.pores('top'))
        geom2['pore.fixed_diameter'] = 1
        Ps = net.pores('top', mode='nor')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps,
                                            throats=net.Ts)
        geom1.add_model(propname='pore.diameter',
                        model=mods.largest_sphere,
                        iters=0)
        assert np.all(np.unique(geom1['pore.diameter']) == [1.5, 5.0])
        geom2['pore.fixed_diameter'] = 6
        geom1.regenerate_models()
        assert np.amin(geom1['pore.diameter']) < 0

    def test_equivalent_diameter(self):
        mod = op.models.geometry.pore_size.equivalent_diameter
        self.geo['pore.volume'] = 1.0
        self.geo.add_model(propname='pore.diameter',
                           model=mod,
                           pore_volume='pore.volume',
                           pore_shape='sphere')
        a = np.unique(self.geo['pore.diameter'])
        b = np.array(1.24070098, ndmin=1)
        assert np.allclose(a, b)
        del self.geo['pore.diameter'], self.geo.models['pore.diameter']
        self.geo.add_model(propname='pore.diameter',
                           model=mod,
                           pore_volume='pore.volume',
                           pore_shape='cube')
        a = np.unique(self.geo['pore.diameter'])
        b = np.array(1.0, ndmin=1)
        assert np.allclose(a, b)


if __name__ == '__main__':

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

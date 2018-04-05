import openpnm as op
import scipy as sp
import scipy.stats as spst
import openpnm.models.geometry.pore_size as mods


class PoreSizeTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.seed'] = sp.rand(self.geo.Np)

    def test_normal(self):
        self.geo.add_model(propname='pore.diameter',
                           model=mods.normal,
                           scale=0.01,
                           loc=0.5,
                           seeds='pore.seed')
        assert 0.45 < sp.mean(self.geo['pore.diameter']) < 0.55
        del self.geo['pore.diameter']

    def test_weibull(self):
        self.geo.add_model(propname='pore.diameter',
                           model=mods.weibull,
                           shape=1.5,
                           scale=0.0001,
                           loc=0.001,
                           seeds='pore.seed')
        assert sp.amin(self.geo['pore.diameter']) > 0.001
        del self.geo['pore.diameter']

    def test_generic(self):
        func = spst.gamma(a=2, loc=0.001, scale=0.0001)
        self.geo.add_model(propname='pore.diameter',
                           model=mods.generic,
                           func=func,
                           seeds='pore.seed')
        assert sp.amin(self.geo['pore.diameter']) > 0.001
        del self.geo['pore.diameter']

    def test_largest_sphere(self):
        net = op.network.Cubic(shape=[5, 5, 5], spacing=[0.1, 0.2, 0.3])
        geo = op.geometry.GenericGeometry(network=net, pores=net.Ps,
                                          throats=net.Ts)
        geo.add_model(propname='pore.diameter',
                      model=mods.largest_sphere,
                      iters=1)
        dmin = sp.amin(geo['pore.diameter'])
        assert dmin <= 0.1
        geo.models['pore.diameter']['iters'] = 5
        geo.regenerate_models()
        assert dmin < sp.amin(geo['pore.diameter'])

    def test_largest_sphere_multiple_geometries(self):
        net = op.network.Cubic(shape=[10, 10, 10], spacing=[5, 5, 5])
        net['pore.coords'][net.pores('top')] += [0, 0, -3]
        geom2 = op.geometry.GenericGeometry(network=net,
                                            pores=net.pores('top'))
        geom2['pore.diameter'] = 1.0
        Ps = net.pores('top', mode='not')
        geom1 = op.geometry.GenericGeometry(network=net, pores=Ps,
                                            throats=net.Ts)
        mod = op.models.geometry.pore_size.largest_sphere
        geom1.add_model(propname='pore.diameter',
                        model=mod,
                        iters=15)
        assert sp.all(geom2['pore.diameter'] == 1.0)
#        assert sp.all(sp.ceil(sp`.unique(geom1['pore.diameter'])) == [3.0, 5.0])


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

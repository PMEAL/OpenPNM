import openpnm as op
import scipy as sp
import openpnm.models.misc as mods
from numpy.testing import assert_approx_equal, assert_array_almost_equal_nulp


class MiscTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)

    def test_constant(self):
        self.geo.add_model(model=mods.constant,
                           propname='pore.value',
                           value=3.3)
        assert sp.all(sp.unique(self.geo['pore.value']) == 3.3)

    def test_product(self):
        self.geo.add_model(model=mods.constant,
                           propname='pore.value1',
                           value=2)
        self.geo.add_model(model=mods.constant,
                           propname='pore.value2',
                           value=2)
        self.geo.add_model(model=mods.product,
                           propname='pore.result1',
                           prop1='pore.value1',
                           prop2='pore.value2')
        assert sp.all(sp.unique(self.geo['pore.result1']) == 4)
        self.geo.add_model(model=mods.constant,
                           propname='pore.value3',
                           value=2)
        self.geo.add_model(model=mods.product,
                           propname='pore.result2',
                           prop1='pore.value1',
                           prop2='pore.value2',
                           prop3='pore.value3')
        assert sp.all(sp.unique(self.geo['pore.result2']) == 8)

    def test_generic_function(self):
        self.geo['pore.rand'] = sp.rand(self.geo.Np)
        self.geo.add_model(model=mods.generic_function,
                           func=sp.clip,
                           propname='pore.clipped',
                           prop='pore.rand',
                           a_min=0.2, a_max=0.8)
        assert sp.amax(self.geo['pore.clipped']) == 0.8
        assert sp.amin(self.geo['pore.clipped']) == 0.2

    def test_scaled(self):
        self.geo['pore.value4'] = 4
        self.geo.add_model(model=mods.scaled,
                           propname='pore.value5',
                           prop='pore.value4',
                           factor=2)
        assert sp.all(sp.unique(self.geo['pore.value5']) == 8)

    def test_linear(self):
        self.geo['pore.value4'] = 4
        self.geo.add_model(model=mods.linear,
                           propname='pore.value6',
                           prop='pore.value4',
                           m=2, b=2)
        assert sp.all(sp.unique(self.geo['pore.value6']) == 10)

    def test_polynomial(self):
        self.geo['pore.value4'] = 4
        self.geo.add_model(model=mods.polynomial,
                           propname='pore.value7',
                           prop='pore.value4',
                           a=[0, 2, 4, 6])
        assert sp.all(sp.unique(self.geo['pore.value7']) == 456)

    def test_random_no_seed(self):
        self.geo.add_model(model=mods.random,
                           propname='pore.seed',
                           element='pore',
                           seed=None)
        temp1 = self.geo['pore.seed'].copy()
        self.geo.regenerate_models()
        temp2 = self.geo['pore.seed'].copy()
        assert sp.all(~(temp1 == temp2))

    def test_random_with_seed(self):
        self.geo.add_model(model=mods.random,
                           propname='pore.seed',
                           element='pore',
                           seed=0)
        temp1 = self.geo['pore.seed'].copy()
        self.geo.regenerate_models()
        temp2 = self.geo['pore.seed'].copy()
        # assert_array_almost_equal_nulp(temp1, temp2)

    def test_random_with_range(self):
        self.geo.add_model(model=mods.random,
                           propname='pore.seed',
                           element='pore',
                           num_range=[0.1, 0.9])
        self.geo.regenerate_models()
        assert sp.amax(self.geo['pore.seed']) <= 0.9
        assert sp.amin(self.geo['pore.seed']) >= 0.1

    def test_from_neighbor_throats_min(self):
        self.geo.pop('pore.seed', None)
        self.geo.models.pop('pore.seed', None)
        self.geo.models.pop('throat.seed', None)
        self.geo['throat.seed'] = sp.rand(self.net.Nt,)
        self.geo.add_model(model=mods.from_neighbor_throats,
                           propname='pore.seed',
                           throat_prop='throat.seed',
                           mode='min')
        assert sp.all(sp.in1d(self.geo['pore.seed'], self.geo['throat.seed']))
        pmax = sp.amax(self.geo['pore.seed'])
        tmax = sp.amax(self.geo['throat.seed'])
        assert pmax <= tmax

    def test_from_neighbor_throats_max(self):
        self.geo.pop('pore.seed', None)
        self.geo.models.pop('pore.seed', None)
        self.geo.models.pop('throat.seed', None)
        self.geo['throat.seed'] = sp.rand(self.net.Nt,)
        self.geo.add_model(model=mods.from_neighbor_throats,
                           propname='pore.seed',
                           throat_prop='throat.seed',
                           mode='max')
        assert sp.all(sp.in1d(self.geo['pore.seed'], self.geo['throat.seed']))
        pmin = sp.amin(self.geo['pore.seed'])
        tmin = sp.amin(self.geo['throat.seed'])
        assert pmin >= tmin

    def test_from_neighbor_throats_mean(self):
        self.geo.pop('pore.seed', None)
        self.geo.models.pop('pore.seed', None)
        self.geo.models.pop('throat.seed', None)
        self.geo['throat.seed'] = sp.rand(self.net.Nt,)
        self.geo.add_model(model=mods.from_neighbor_throats,
                           propname='pore.seed',
                           throat_prop='throat.seed',
                           mode='mean')
        tmax = sp.amax(self.geo['throat.seed'])
        tmin = sp.amin(self.geo['throat.seed'])
        assert sp.all(self.geo['pore.seed'] > tmin)
        assert sp.all(self.geo['pore.seed'] < tmax)

    def test_from_neighbor_pores_min(self):
        self.geo.remove_model('throat.seed')
        self.geo['pore.seed'] = sp.rand(self.net.Np,)
        self.geo.add_model(model=mods.from_neighbor_pores,
                           propname='throat.seed',
                           pore_prop='pore.seed',
                           mode='min')
        P12 = self.net['throat.conns']
        tseed = sp.amin(self.geo['pore.seed'][P12], axis=1)
        assert_array_almost_equal_nulp(self.geo['throat.seed'], tseed)

    def test_from_neighbor_pores_max(self):
        self.geo.remove_model('throat.seed')
        self.geo['pore.seed'] = sp.rand(self.net.Np,)
        self.geo.add_model(model=mods.from_neighbor_pores,
                           propname='throat.seed',
                           pore_prop='pore.seed',
                           mode='max')
        P12 = self.net['throat.conns']
        tseed = sp.amax(self.geo['pore.seed'][P12], axis=1)
        assert_array_almost_equal_nulp(self.geo['throat.seed'], tseed)

    def test_from_neighbor_pores_mean(self):
        self.geo.remove_model('throat.seed')
        self.geo['pore.seed'] = sp.rand(self.net.Np,)
        self.geo.add_model(model=mods.from_neighbor_pores,
                           propname='throat.seed',
                           pore_prop='pore.seed',
                           mode='mean')
        P12 = self.net['throat.conns']
        tseed = sp.mean(self.geo['pore.seed'][P12], axis=1)
        assert_array_almost_equal_nulp(self.geo['throat.seed'], tseed)

    def test_from_neighbors_multi_geom(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        Ps1 = net.pores('internal')
        Ts1 = net.throats('internal')
        geo1 = op.geometry.GenericGeometry(network=net,
                                           pores=Ps1,
                                           throats=Ts1)
        Ps2 = net.pores('internal', mode='not')
        Ts2 = net.throats('internal', mode='not')
        geo2 = op.geometry.GenericGeometry(network=net,
                                           pores=Ps2,
                                           throats=Ts2)
        geo1['pore.rand1'] = sp.random.random(geo1.Np)
        geo2['pore.rand1'] = sp.random.random(geo2.Np)
        geo1.add_model(model=mods.from_neighbor_pores,
                       propname='throat.rand1',
                       pore_prop='pore.rand1',
                       mode='min')
        test = sp.amin(net['pore.rand1'][net['throat.conns']], axis=1)[Ts1]
        assert sp.all(test == geo1['throat.rand1'])
        geo1['throat.rand2'] = sp.random.random(geo1.Nt)
        geo2['throat.rand2'] = sp.random.random(geo2.Nt)
        geo2.add_model(model=mods.from_neighbor_throats,
                       propname='pore.rand2',
                       throat_prop='throat.rand2',
                       mode='max')
        test = sp.zeros(geo2.Np).astype(bool)
        for i, pore in enumerate(net.pores(geo2.name)):
            Ts = net.find_neighbor_throats(pores=pore)
            T_max = sp.amax(net['throat.rand2'][Ts])
            test[i] = net['pore.rand2'][pore] == T_max
        assert sp.all(test)

if __name__ == '__main__':

    t = MiscTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

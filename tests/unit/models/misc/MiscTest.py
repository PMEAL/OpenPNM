import numpy as np
import openpnm as op
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
        assert np.all(np.unique(self.geo['pore.value']) == 3.3)

    def test_product(self):
        self.geo.add_model(model=mods.constant,
                           propname='pore.value1',
                           value=2)
        self.geo.add_model(model=mods.constant,
                           propname='pore.value2',
                           value=2)
        self.geo.add_model(model=mods.product,
                           propname='pore.result1',
                           props=['pore.value1', 'pore.value2'])
        assert np.all(np.unique(self.geo['pore.result1']) == 4)
        self.geo.add_model(model=mods.constant,
                           propname='pore.value3',
                           value=2)
        self.geo.add_model(model=mods.product,
                           propname='pore.result2',
                           props=['pore.value1', 'pore.value2', 'pore.value3'])
        assert np.all(np.unique(self.geo['pore.result2']) == 8)

    def test_generic_function(self):
        self.geo['pore.rand'] = np.random.rand(self.geo.Np)
        self.geo.add_model(model=mods.generic_function,
                           func=np.clip,
                           propname='pore.clipped',
                           prop='pore.rand',
                           a_min=0.2, a_max=0.8)
        assert np.amax(self.geo['pore.clipped']) == 0.8
        assert np.amin(self.geo['pore.clipped']) == 0.2

    def test_scaled(self):
        self.geo['pore.value4'] = 4
        self.geo.add_model(model=mods.scaled,
                           propname='pore.value5',
                           prop='pore.value4',
                           factor=2)
        assert np.all(np.unique(self.geo['pore.value5']) == 8)

    def test_linear(self):
        self.geo['pore.value4'] = 4
        self.geo.add_model(model=mods.linear,
                           propname='pore.value6',
                           prop='pore.value4',
                           m=2, b=2)
        assert np.all(np.unique(self.geo['pore.value6']) == 10)

    def test_polynomial(self):
        self.geo['pore.value4'] = 4
        self.geo.add_model(model=mods.polynomial,
                           propname='pore.value7',
                           prop='pore.value4',
                           a=[0, 2, 4, 6])
        assert np.all(np.unique(self.geo['pore.value7']) == 456)

    def test_random_no_seed(self):
        self.geo.add_model(model=mods.random,
                           propname='pore.seed',
                           element='pore',
                           seed=None)
        temp1 = self.geo['pore.seed'].copy()
        self.geo.regenerate_models()
        temp2 = self.geo['pore.seed'].copy()
        assert np.all(~(temp1 == temp2))

    def test_random_with_seed(self):
        self.geo.add_model(model=mods.random,
                           propname='pore.seed',
                           element='pore',
                           seed=0)
        temp1 = self.geo['pore.seed'].copy()
        self.geo.regenerate_models()
        temp2 = self.geo['pore.seed'].copy()
        assert_array_almost_equal_nulp(temp1, temp2)

    def test_random_with_range(self):
        self.geo.add_model(model=mods.random,
                           propname='pore.seed',
                           element='pore',
                           num_range=[0.1, 0.9])
        self.geo.regenerate_models()
        assert np.amax(self.geo['pore.seed']) <= 0.9
        assert np.amin(self.geo['pore.seed']) >= 0.1

    def test_from_neighbor_throats_min(self):
        self.geo.pop('pore.seed', None)
        self.geo.models.pop('pore.seed', None)
        self.geo.models.pop('throat.seed', None)
        self.geo['throat.seed'] = np.linspace(0, 1, self.net.Nt)
        self.geo.add_model(model=mods.from_neighbor_throats,
                           propname='pore.seed',
                           prop='throat.seed',
                           mode='min')
        assert np.all(np.in1d(self.geo['pore.seed'], self.geo['throat.seed']))
        assert np.isclose(self.geo['throat.seed'].mean(), 0.5)
        assert np.isclose(self.geo['pore.seed'].mean(), 0.16454849498327762)

    def test_from_neighbor_throats_max(self):
        self.geo.pop('pore.seed', None)
        self.geo.models.pop('pore.seed', None)
        self.geo.models.pop('throat.seed', None)
        self.geo['throat.seed'] = np.linspace(0, 1, self.net.Nt)
        self.geo.add_model(model=mods.from_neighbor_throats,
                           propname='pore.seed',
                           prop='throat.seed',
                           mode='max')
        assert np.all(np.in1d(self.geo['pore.seed'], self.geo['throat.seed']))
        assert np.isclose(self.geo['throat.seed'].mean(), 0.5)
        assert np.isclose(self.geo['pore.seed'].mean(), 0.8595317725752508)

    def test_from_neighbor_throats_mean(self):
        self.geo.pop('pore.seed', None)
        self.geo.models.pop('pore.seed', None)
        self.geo.models.pop('throat.seed', None)
        self.geo['throat.seed'] = np.linspace(0, 1, self.net.Nt)
        self.geo.add_model(model=mods.from_neighbor_throats,
                           propname='pore.seed',
                           prop='throat.seed',
                           mode='mean')
        assert np.isclose(self.geo['throat.seed'].mean(), 0.5)
        assert np.isclose(self.geo['pore.seed'].mean(), 0.5)

    def test_neighbor_pores_with_nans(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        net['pore.values'] = 1.0
        net['pore.values'][0] = np.nan
        f = mods.from_neighbor_pores
        with_nans = f(target=net, prop='pore.values',
                      ignore_nans=False, mode='min')
        assert np.any(np.isnan(with_nans))
        no_nans = f(target=net, prop='pore.values',
                    ignore_nans=True, mode='min')
        assert np.all(~np.isnan(no_nans))
        with_nans = f(target=net, prop='pore.values',
                      ignore_nans=False, mode='max')
        assert np.any(np.isnan(with_nans))
        no_nans = f(target=net, prop='pore.values',
                    ignore_nans=True, mode='max')
        assert np.all(~np.isnan(no_nans))
        with_nans = f(target=net, prop='pore.values',
                      ignore_nans=False, mode='mean')
        assert np.any(np.isnan(with_nans))
        no_nans = f(target=net, prop='pore.values',
                    ignore_nans=True, mode='mean')
        assert np.all(~np.isnan(no_nans))

    def test_neighbor_throats_mode_min_with_nans(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        net['throat.values'] = np.linspace(0, 1, net.Nt)
        net['throat.values'][0] = np.nan
        f = mods.from_neighbor_throats
        with_nans = f(target=net, prop='throat.values',
                      ignore_nans=False, mode='min')
        assert np.any(np.isnan(with_nans))
        no_nans = f(target=net, prop='throat.values',
                    ignore_nans=True, mode='min')
        assert np.all(~np.isnan(no_nans))
        assert np.all(~np.isinf(no_nans))
        assert np.allclose(no_nans, np.array([0.36363636, 0.45454545,
                                              0.09090909, 0.09090909,
                                              0.18181818, 0.18181818,
                                              0.27272727, 0.27272727]))

    def test_neighbor_throats_mode_max_with_nans(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        net['throat.values'] = np.linspace(0, 1, net.Nt)
        net['throat.values'][0] = np.nan
        f = mods.from_neighbor_throats
        with_nans = f(target=net, prop='throat.values',
                      ignore_nans=False, mode='max')
        assert np.any(np.isnan(with_nans))
        no_nans = f(target=net, prop='throat.values',
                    ignore_nans=True, mode='max')
        assert np.all(~np.isnan(no_nans))
        assert np.all(~np.isinf(no_nans))
        assert np.allclose(no_nans, np.array([0.72727273, 0.81818182,
                                              0.90909091, 1.00000000,
                                              0.72727273, 0.81818182,
                                              0.90909091, 1.00000000]))

    def test_neighbor_throats_mode_mean_with_nans(self):
        net = op.network.Cubic(shape=[2, 2, 2])
        net['throat.values'] = np.linspace(0, 1, net.Nt)
        net['throat.values'][0] = np.nan
        f = mods.from_neighbor_throats
        with_nans = f(target=net, prop='throat.values',
                      ignore_nans=False, mode='mean')
        assert np.any(np.isnan(with_nans))
        no_nans = f(target=net, prop='throat.values',
                    ignore_nans=True, mode='mean')
        assert np.all(~np.isnan(no_nans))
        assert np.all(~np.isinf(no_nans))
        assert np.allclose(no_nans, np.array([0.54545455, 0.63636364,
                                              0.45454545, 0.51515152,
                                              0.48484848, 0.54545455,
                                              0.57575758, 0.63636364]))

    def test_from_neighbor_pores_min(self):
        self.geo.remove_model('throat.seed')
        self.geo['pore.seed'] = np.random.rand(self.net.Np,)
        self.geo.add_model(model=mods.from_neighbor_pores,
                           propname='throat.seed',
                           prop='pore.seed',
                           mode='min')
        P12 = self.net['throat.conns']
        tseed = np.amin(self.geo['pore.seed'][P12], axis=1)
        assert_array_almost_equal_nulp(self.geo['throat.seed'], tseed)

    def test_from_neighbor_pores_max(self):
        self.geo.remove_model('throat.seed')
        self.geo['pore.seed'] = np.random.rand(self.net.Np,)
        self.geo.add_model(model=mods.from_neighbor_pores,
                           propname='throat.seed',
                           prop='pore.seed',
                           mode='max')
        P12 = self.net['throat.conns']
        tseed = np.amax(self.geo['pore.seed'][P12], axis=1)
        assert_array_almost_equal_nulp(self.geo['throat.seed'], tseed)

    def test_from_neighbor_pores_mean(self):
        self.geo.remove_model('throat.seed')
        self.geo['pore.seed'] = np.random.rand(self.net.Np,)
        self.geo.add_model(model=mods.from_neighbor_pores,
                           propname='throat.seed',
                           prop='pore.seed',
                           mode='mean')
        P12 = self.net['throat.conns']
        tseed = np.mean(self.geo['pore.seed'][P12], axis=1)
        assert_array_almost_equal_nulp(self.geo['throat.seed'], tseed)

    def test_from_neighbors_multi_geom(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        net.add_boundary_pores()
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
        geo1['pore.rand1'] = np.random.random(geo1.Np)
        geo2['pore.rand1'] = np.random.random(geo2.Np)
        geo1.add_model(model=mods.from_neighbor_pores,
                       propname='throat.rand1',
                       prop='pore.rand1',
                       mode='min')
        test = np.amin(net['pore.rand1'][net['throat.conns']], axis=1)[Ts1]
        assert np.all(test == geo1['throat.rand1'])
        geo1['throat.rand2'] = np.random.random(geo1.Nt)
        geo2['throat.rand2'] = np.random.random(geo2.Nt)
        geo2.add_model(model=mods.from_neighbor_throats,
                       propname='pore.rand2',
                       prop='throat.rand2',
                       mode='max')
        test = np.zeros(geo2.Np).astype(bool)
        for i, pore in enumerate(net.pores(geo2.name)):
            Ts = net.find_neighbor_throats(pores=pore)
            T_max = np.amax(net['throat.rand2'][Ts])
            test[i] = net['pore.rand2'][pore] == T_max
        assert np.all(test)

    def test_invert(self):
        net = op.network.Cubic(shape=[5, 5, 5])
        net['pore.diameter'] = 2.0
        net.add_model(propname='pore.entry_pressure',
                      prop='pore.diameter',
                      model=mods.basic_math.invert)
        assert net['pore.entry_pressure'][0] == 0.5


if __name__ == '__main__':

    t = MiscTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

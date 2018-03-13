import openpnm as op
import scipy as sp
import openpnm.models.misc as mods


class PoreMiscTest:

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
                           value=4)
        self.geo.add_model(model=mods.product,
                           propname='pore.value3',
                           props=['pore.value1', 'pore.value2'])
        assert sp.all(sp.unique(self.geo['pore.value3']) == 8)

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
#        assert np.testing.assert_array_almost_equal_nulp(temp1, temp2)

    def test_random_with_range(self):
        self.geo.add_model(model=mods.random,
                           propname='pore.seed',
                           element='pore',
                           num_range=[0.1, 0.9])
        self.geo.regenerate_models()
        assert sp.amax(self.geo['pore.seed']) <= 0.9
        assert sp.amin(self.geo['pore.seed']) >= 0.1

    def test_from_neighbor_throats_min(self):
        catch = self.geo.pop('pore.seed', None)
        catch = self.geo.models.pop('pore.seed', None)
        catch = self.geo.models.pop('throat.seed', None)
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
        catch = self.geo.pop('pore.seed', None)
        catch = self.geo.models.pop('pore.seed', None)
        catch = self.geo.models.pop('throat.seed', None)
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
        catch = self.geo.pop('pore.seed', None)
        catch = self.geo.models.pop('pore.seed', None)
        catch = self.geo.models.pop('throat.seed', None)
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
        assert sp.allclose(self.geo['throat.seed'], tseed)

    def test_from_neighbor_pores_max(self):
        self.geo.remove_model('throat.seed')
        self.geo['pore.seed'] = sp.rand(self.net.Np,)
        self.geo.add_model(model=mods.from_neighbor_pores,
                           propname='throat.seed',
                           pore_prop='pore.seed',
                           mode='max')
        P12 = self.net['throat.conns']
        tseed = sp.amax(self.geo['pore.seed'][P12], axis=1)
        assert sp.allclose(self.geo['throat.seed'], tseed)

    def test_from_neighbor_pores_mean(self):
        self.geo.remove_model('throat.seed')
        self.geo['pore.seed'] = sp.rand(self.net.Np,)
        self.geo.add_model(model=mods.from_neighbor_pores,
                           propname='throat.seed',
                           pore_prop='pore.seed',
                           mode='mean')
        P12 = self.net['throat.conns']
        tseed = sp.mean(self.geo['pore.seed'][P12], axis=1)
        assert sp.allclose(self.geo['throat.seed'], tseed)


if __name__ == '__main__':

    t = PoreMiscTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

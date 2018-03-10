import openpnm as op
import scipy as sp
import openpnm.models.geometry as gm


class PoreMiscTest:
    def test_constant(self):
        pass

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)

    def test_random_no_seed(self):
        mod = gm.pore_misc.random
        self.geo.add_model(model=mod,
                           propname='pore.seed',
                           seed=None)
        temp1 = self.geo['pore.seed'].copy()
        self.geo.regenerate_models()
        temp2 = self.geo['pore.seed'].copy()
        assert sp.all(~(temp1 == temp2))

    def test_random_with_seed(self):
        mod = gm.pore_misc.random
        self.geo.add_model(model=mod,
                           propname='pore.seed',
                           seed=0)
        temp1 = self.geo['pore.seed'].copy()
        self.geo.regenerate_models()
        temp2 = self.geo['pore.seed'].copy()
#        assert np.testing.assert_array_almost_equal_nulp(temp1, temp2)

    def test_random_with_range(self):
        mod = gm.pore_misc.random
        self.geo.add_model(model=mod,
                           propname='pore.seed',
                           num_range=[0.1, 0.9])
        self.geo.regenerate_models()
        assert sp.amax(self.geo['pore.seed']) <= 0.9
        assert sp.amin(self.geo['pore.seed']) >= 0.1

    def test_neighbor_min(self):
        catch = self.geo.pop('pore.seed', None)
        catch = self.geo.models.pop('pore.seed', None)
        catch = self.geo.models.pop('throat.seed', None)
        mod = gm.pore_misc.neighbor
        self.geo['throat.seed'] = sp.rand(self.net.Nt,)
        self.geo.add_model(model=mod,
                           propname='pore.seed',
                           throat_prop='throat.seed',
                           mode='min')
        assert sp.all(sp.in1d(self.geo['pore.seed'], self.geo['throat.seed']))
        pmax = sp.amax(self.geo['pore.seed'])
        tmax = sp.amax(self.geo['throat.seed'])
        assert pmax <= tmax

    def test_neighbor_max(self):
        catch = self.geo.pop('pore.seed', None)
        catch = self.geo.models.pop('pore.seed', None)
        catch = self.geo.models.pop('throat.seed', None)
        mod = gm.pore_misc.neighbor
        self.geo['throat.seed'] = sp.rand(self.net.Nt,)
        self.geo.add_model(model=mod,
                           propname='pore.seed',
                           throat_prop='throat.seed',
                           mode='max')
        assert sp.all(sp.in1d(self.geo['pore.seed'], self.geo['throat.seed']))
        pmin = sp.amin(self.geo['pore.seed'])
        tmin = sp.amin(self.geo['throat.seed'])
        assert pmin >= tmin

    def test_neighbor_mean(self):
        catch = self.geo.pop('pore.seed', None)
        catch = self.geo.models.pop('pore.seed', None)
        catch = self.geo.models.pop('throat.seed', None)
        mod = gm.pore_misc.neighbor
        self.geo['throat.seed'] = sp.rand(self.net.Nt,)
        self.geo.add_model(model=mod,
                           propname='pore.seed',
                           throat_prop='throat.seed',
                           mode='mean')
        tmax = sp.amax(self.geo['throat.seed'])
        tmin = sp.amin(self.geo['throat.seed'])
        assert sp.all(self.geo['pore.seed'] > tmin)
        assert sp.all(self.geo['pore.seed'] < tmax)


if __name__ == '__main__':

    t = PoreMiscTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

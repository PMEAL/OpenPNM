import OpenPNM
import scipy as sp
import OpenPNM.Geometry.models as gm


class ThroatMiscTest:

    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts)

    def test_random_no_seed(self):
        mod = gm.throat_misc.random
        self.geo.models.add(model=mod,
                            propname='throat.seed',
                            seed=None)
        temp1 = self.geo['throat.seed'].copy()
        self.geo.models.regenerate()
        temp2 = self.geo['throat.seed'].copy()
        assert sp.all(~(temp1 == temp2))

    def test_random_with_seed(self):
        mod = gm.throat_misc.random
        self.geo.models.add(model=mod,
                            propname='throat.seed',
                            seed=0)
        temp1 = self.geo['throat.seed'].copy()
        self.geo.models.regenerate()
        temp2 = self.geo['throat.seed'].copy()
        assert sp.all(temp1 == temp2)

    def test_random_with_range(self):
        mod = gm.throat_misc.random
        self.geo.models.add(model=mod,
                            propname='throat.seed',
                            num_range=[0.1, 0.9])
        assert sp.amax(self.geo['throat.seed']) <= 0.9
        assert sp.amin(self.geo['throat.seed']) >= 0.1

    def test_neighbor_min(self):
        catch = self.geo.pop('throat.seed', None)
        catch = self.geo.models.pop('throat.seed', None)
        mod = gm.throat_misc.neighbor
        self.geo['pore.seed'] = sp.rand(self.net.Np,)
        self.geo.models.add(model=mod,
                            propname='throat.seed',
                            pore_prop='pore.seed',
                            mode='min')
        P12 = self.net['throat.conns']
        tseed = sp.amin(self.geo['pore.seed'][P12], axis=1)
        assert sp.allclose(self.geo['throat.seed'], tseed)

    def test_neighbor_max(self):
        catch = self.geo.pop('throat.seed', None)
        catch = self.geo.models.pop('throat.seed', None)
        mod = gm.throat_misc.neighbor
        self.geo['pore.seed'] = sp.rand(self.net.Np,)
        self.geo.models.add(model=mod,
                            propname='throat.seed',
                            pore_prop='pore.seed',
                            mode='max')
        P12 = self.net['throat.conns']
        tseed = sp.amax(self.geo['pore.seed'][P12], axis=1)
        assert sp.allclose(self.geo['throat.seed'], tseed)

    def test_neighbor_mean(self):
        catch = self.geo.pop('throat.seed', None)
        catch = self.geo.models.pop('throat.seed', None)
        mod = gm.throat_misc.neighbor
        self.geo['pore.seed'] = sp.rand(self.net.Np,)
        self.geo.models.add(model=mod,
                            propname='throat.seed',
                            pore_prop='pore.seed',
                            mode='mean')
        P12 = self.net['throat.conns']
        tseed = sp.mean(self.geo['pore.seed'][P12], axis=1)
        assert sp.allclose(self.geo['throat.seed'], tseed)

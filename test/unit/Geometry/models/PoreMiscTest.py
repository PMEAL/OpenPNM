import OpenPNM
import scipy as sp
import OpenPNM.Geometry.models as gm


class PoreMiscTest:
    def test_constant(self):
        pass

    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts)

    def test_random_no_seed(self):
        mod = gm.pore_misc.random
        self.geo.models.add(model=mod,
                            propname='pore.seed',
                            seed=None)
        temp1 = self.geo['pore.seed'].copy()
        self.geo.models.regenerate()
        temp2 = self.geo['pore.seed'].copy()
        assert sp.all(~(temp1 == temp2))

    def test_random_with_seed(self):
        mod = gm.pore_misc.random
        self.geo.models.add(model=mod,
                            propname='pore.seed',
                            seed=0)
        temp1 = self.geo['pore.seed'].copy()
        self.geo.models.regenerate()
        temp2 = self.geo['pore.seed'].copy()
        assert sp.all(temp1 == temp2)

    def test_random_with_range(self):
        mod = gm.pore_misc.random
        self.geo.models.add(model=mod,
                            propname='pore.seed',
                            num_range=[0.1, 0.9])
        assert sp.amax(self.geo['pore.seed']) <= 0.9
        assert sp.amin(self.geo['pore.seed']) >= 0.1

    def test_neighbor_min(self):
        catch = self.geo.pop('pore.seed', None)
        catch = self.geo.models.pop('pore.seed', None)
        catch = self.geo.models.pop('throat.seed', None)
        mod = gm.pore_misc.neighbor
        self.geo['throat.seed'] = sp.rand(self.net.Nt,)
        self.geo.models.add(model=mod,
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
        self.geo.models.add(model=mod,
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
        self.geo.models.add(model=mod,
                            propname='pore.seed',
                            throat_prop='throat.seed',
                            mode='mean')
        tmax = sp.amax(self.geo['throat.seed'])
        tmin = sp.amin(self.geo['throat.seed'])
        assert sp.all(self.geo['pore.seed'] > tmin)
        assert sp.all(self.geo['pore.seed'] < tmax)

import OpenPNM
import scipy as sp


class ThroatSeedTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts)
        self.geo['pore.seed'] = sp.rand(self.net.Np)

    def test_random(self):
        self.geo.models.add(propname='throat.seed',
                            model=OpenPNM.Geometry.models.throat_seed.random,
                            seed=0,
                            num_range=[0.1, 2])
        assert sp.amax(self.geo['throat.seed']) > 1.9
        assert sp.amin(self.geo['throat.seed']) > 0.1

    def test_neighbor(self):
        self.geo.models.add(propname='throat.seed_max',
                            model=OpenPNM.Geometry.models.throat_seed.neighbor,
                            mode='max')
        self.geo.models.add(propname='throat.seed_min',
                            model=OpenPNM.Geometry.models.throat_seed.neighbor,
                            mode='min')
        self.geo.models.add(propname='throat.seed_mean',
                            model=OpenPNM.Geometry.models.throat_seed.neighbor,
                            mode='mean')
        assert sp.all(self.geo['throat.seed_min'] <= self.geo['throat.seed_max'])
        assert sp.all(self.geo['throat.seed_min'] <= self.geo['throat.seed_mean'])
        assert sp.all(self.geo['throat.seed_mean'] <= self.geo['throat.seed_max'])

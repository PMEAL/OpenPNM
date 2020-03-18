import numpy as np
import scipy as sp
import openpnm as op
import openpnm.models.geometry.throat_seed as mods


class ThroatSeedTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.seed'] = np.random.rand(self.net.Np)

    def test_random(self):
        self.geo.add_model(propname='throat.seed',
                           model=mods.random,
                           seed=0,
                           num_range=[0.1, 2])
        self.geo.regenerate_models()
        assert np.amax(self.geo['throat.seed']) > 1.9
        assert np.amin(self.geo['throat.seed']) > 0.1

    def test_neighbor(self):
        self.geo.add_model(propname='throat.seed_max',
                           model=mods.from_neighbor_pores,
                           mode='max')
        self.geo.add_model(propname='throat.seed_min',
                           model=mods.from_neighbor_pores,
                           mode='min')
        self.geo.add_model(propname='throat.seed_mean',
                           model=mods.from_neighbor_pores,
                           mode='mean')
        self.geo.regenerate_models()
        assert np.all(self.geo['throat.seed_min'] <= self.geo['throat.seed_max'])
        assert np.all(self.geo['throat.seed_min'] <= self.geo['throat.seed_mean'])
        assert np.all(self.geo['throat.seed_mean'] <= self.geo['throat.seed_max'])


if __name__ == '__main__':

    t = ThroatSeedTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

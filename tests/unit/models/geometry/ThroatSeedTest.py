import numpy as np
import openpnm as op
import openpnm.models.geometry.throat_seed as mods


class ThroatSeedTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.net['pore.seed'] = np.random.rand(self.net.Np)

    def test_random(self):
        self.net.add_model(propname='throat.seed',
                           model=mods.random,
                           seed=0,
                           num_range=[0.1, 2])
        self.net.regenerate_models()
        assert np.amax(self.net['throat.seed']) > 1.9
        assert np.amin(self.net['throat.seed']) > 0.1

    def test_neighbor(self):
        self.net.add_model(propname='throat.seed_max',
                           model=mods.from_neighbor_pores,
                           mode='max')
        self.net.add_model(propname='throat.seed_min',
                           model=mods.from_neighbor_pores,
                           mode='min')
        self.net.add_model(propname='throat.seed_mean',
                           model=mods.from_neighbor_pores,
                           mode='mean')
        self.net.regenerate_models()
        assert np.all(self.net['throat.seed_min'] <= self.net['throat.seed_max'])
        assert np.all(self.net['throat.seed_min'] <= self.net['throat.seed_mean'])
        assert np.all(self.net['throat.seed_mean'] <= self.net['throat.seed_max'])


if __name__ == '__main__':

    t = ThroatSeedTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

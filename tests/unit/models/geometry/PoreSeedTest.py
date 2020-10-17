import numpy as np
import scipy as sp
import openpnm as op
import openpnm.models.geometry.pore_seed as mods


class PoreSeedTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)

    def test_random(self):
        f = mods.random
        self.geo.add_model(propname='pore.seed',
                           model=f,
                           seed=0,
                           num_range=[0.1, 2])
        assert np.amax(self.geo['pore.seed']) > 1
        assert np.amin(self.geo['pore.seed']) < 1

    def test_spatially_correlated(self):
        f = mods.spatially_correlated
        self.geo.add_model(propname='pore.seed',
                           model=f,
                           weights=[2, 2, 2],
                           regen_mode='normal')
        assert np.amin(self.geo['pore.seed'] > 0)
        assert np.amax(self.geo['pore.seed'] < 1)

    def test_spatially_correlated_zero_weights(self):
        f = mods.spatially_correlated
        self.geo.add_model(propname='pore.seed',
                           model=f,
                           weights=[0, 0, 0],
                           regen_mode='normal')
        assert np.amin(self.geo['pore.seed'] > 0)
        assert np.amax(self.geo['pore.seed'] < 1)


if __name__ == '__main__':

    t = PoreSeedTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

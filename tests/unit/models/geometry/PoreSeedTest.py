import openpnm as op
import scipy as sp
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
        assert sp.amax(self.geo['pore.seed']) > 1
        assert sp.amin(self.geo['pore.seed']) < 1

    def test_spatially_correlated(self):
        f = mods.spatially_correlated
        self.geo.add_model(propname='pore.seed',
                           model=f,
                           weights=[2, 2, 2],
                           regen_mode='normal')
        assert sp.amin(self.geo['pore.seed'] > 0)
        assert sp.amax(self.geo['pore.seed'] < 1)

    def test_spatially_correlated_zero_weights(self):
        f = mods.spatially_correlated
        self.geo.add_model(propname='pore.seed',
                           model=f,
                           weights=[0, 0, 0],
                           regen_mode='normal')
        assert sp.amin(self.geo['pore.seed'] > 0)
        assert sp.amax(self.geo['pore.seed'] < 1)

    def test_location_adjusted(self):
        image1 = sp.ones([5, 5, 5])
        image1[:, :, 0] = 0.5
        f = mods.location_adjusted
        self.geo.add_model(propname='pore.seed',
                           model=f,
                           image=image1)
        Ps = self.net.pores('bottom')
        a = sp.mean(self.geo["pore.seed"][Ps])
        Other = self.net.pores('bottom', mode='not')
        b = sp.mean(self.geo["pore.seed"][Other])
        assert a < b
        image2 = sp.ones([2, 2])
        self.geo.add_model(propname='pore.seed',
                           model=f,
                           image=image2,
                           regen_mode='normal')
        self.geo.add_model(propname='pore.seed',
                           model=f,
                           image=-image1,
                           regen_mode='normal')


if __name__ == '__main__':

    t = PoreSeedTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

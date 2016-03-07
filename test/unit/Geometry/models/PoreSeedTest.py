import OpenPNM
import scipy as sp


class PoreSeedTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts)

    def test_random(self):
        f = OpenPNM.Geometry.models.pore_seed.random
        self.geo.models.add(propname='pore.seed',
                            model=f,
                            seed=0,
                            num_range=[0.1, 2])
        assert sp.amax(self.geo['pore.seed']) > 1.9
        assert sp.amin(self.geo['pore.seed']) > 0.1

    def test_spatially_correlated(self):
        f = OpenPNM.Geometry.models.pore_seed.spatially_correlated
        self.geo.models.add(propname='pore.seed',
                            model=f,
                            weights=[2, 2, 2])
        assert sp.amin(self.geo['pore.seed'] > 0)
        assert sp.amax(self.geo['pore.seed'] < 1)

    def test_spatially_correlated_zero_weights(self):
        f = OpenPNM.Geometry.models.pore_seed.spatially_correlated
        self.geo.models.add(propname='pore.seed',
                            model=f,
                            weights=[0, 0, 0])
        assert sp.amin(self.geo['pore.seed'] > 0)
        assert sp.amax(self.geo['pore.seed'] < 1)

    def test_location_adjusted(self):
        image1 = sp.ones([5, 5, 5])
        image1[:, :, 0] = 0.5
        f = OpenPNM.Geometry.models.pore_seed.location_adjusted
        self.geo.models.add(propname='pore.seed',
                            model=f,
                            image=image1)
        Ps = self.net.pores('bottom')
        a = sp.mean(self.geo["pore.seed"][Ps])
        Other = self.net.pores('bottom', mode='not')
        b = sp.mean(self.geo["pore.seed"][Other])
        assert a < b
        image2 = sp.ones([2, 2])
        self.geo.models.add(propname='pore.seed',
                            model=f,
                            image=image2)
        self.geo.models.add(propname='pore.seed',
                            model=f,
                            image=-image1)

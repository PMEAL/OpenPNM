import openpnm as op
import openpnm.models.geometry.throat_length as mods
import numpy as np
from numpy.testing import assert_approx_equal


class ThroatStraightTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.air,
                                              geometry=self.geo)
        self.geo['pore.diameter'] = .05

    def test_straight(self):
        self.geo.add_model(propname='throat.length',
                           model=mods.straight,
                           regen_mode='normal')
        a = np.array([0.95])
        b = np.unique(self.geo['throat.length'])
        assert np.allclose(a, b)
        self.geo['pore.diameter'] = 1.1
        self.geo.regenerate_models(propnames=['throat.length'])
        assert np.amin(self.geo['throat.length']) == 1e-9

    def test_spherical_pores(self):
        self.geo['throat.diameter'] = 0.035
        self.geo['pore.diameter'] = 0.05
        self.geo.add_model(propname='throat.conduit_lengths',
                           model=mods.spherical_pores,
                           regen_mode='normal')
        L1 = np.unique(self.geo['throat.conduit_lengths.pore1'])
        Lt = np.unique(self.geo['throat.conduit_lengths.throat'])
        L2 = np.unique(self.geo['throat.conduit_lengths.pore2'])
        assert np.allclose(0.01785357, L1)
        assert np.allclose(0.96429286, Lt)
        assert np.allclose(0.01785357, L2)

    def test_circular_pores(self):
        self.geo['throat.diameter'] = 0.035
        self.geo['pore.diameter'] = 0.05
        self.geo.add_model(propname='throat.conduit_lengths',
                           model=mods.circular_pores,
                           regen_mode='normal')
        L1 = np.unique(self.geo['throat.conduit_lengths.pore1'])
        Lt = np.unique(self.geo['throat.conduit_lengths.throat'])
        L2 = np.unique(self.geo['throat.conduit_lengths.pore2'])
        assert np.allclose(0.01785357, L1)
        assert np.allclose(0.96429286, Lt)
        assert np.allclose(0.01785357, L2)

if __name__ == '__main__':

    t = ThroatStraightTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

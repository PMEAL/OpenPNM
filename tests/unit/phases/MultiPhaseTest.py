import openpnm as op
import scipy as sp
import pytest


class MultiPhaseTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.water = op.phases.Water(network=self.net)
        self.air = op.phases.Air(network=self.net)

    def test_multiphase_init(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.air,
                                                           self.water])
        assert sp.all(m['pore.occupancy.all'] == 0.0)
        assert sp.all(m['throat.occupancy.all'] == 0.0)
        assert self.air.name in m.settings['phases']
        assert self.water.name in m.settings['phases']

    def test_multiphase_no_occupancy_yet(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.air,
                                                           self.water])
        self.water['pore.temperature'] = 300
        assert sp.all(self.water['pore.temperature'] == 300)
        with pytest.raises(Exception):
            m['pore.temperature']

    def test_multiphase_occupancy_set_single_phase(self):
        m = op.phases.MultiPhase(network=self.net)
        self.water['pore.temperature'] = 300
        assert sp.all(self.water['pore.temperature'] == 300)
        m.set_occupancy(phase=self.water, Pvals=1, Tvals=1)
        assert sp.all(m['pore.temperature'] == 300)

    def test_multiphase_occupancy_set_two_phase(self):
        m = op.phases.MultiPhase(network=self.net)
        self.water['pore.temperature'] = 300
        self.air['pore.temperature'] = 200
        assert sp.all(self.water['pore.temperature'] == 300)
        assert sp.all(self.air['pore.temperature'] == 200)
        Ps = self.net['pore.coords'][:, 0] < 3
        Ts = self.net.tomask(throats=self.net.find_neighbor_throats(Ps))
        m.set_occupancy(phase=self.water, Pvals=Ps, Tvals=Ts)
        m.set_occupancy(phase=self.air, Pvals=~Ps, Tvals=~Ts)
        assert sp.all(m['pore.temperature'] >= 200)
        assert sp.all(m['pore.temperature'] <= 300)

    def test_multiphase_invalid_phase(self):
        pn = op.network.Cubic(shape=[3, 3, 3])
        water = op.phases.Water(network=pn)
        m = op.phases.MultiPhase(network=self.net)
        with pytest.raises(Exception):
            m.set_occupancy(phase=water)

    def test_multiphase_invalid_occupancy(self):
        m = op.phases.MultiPhase(network=self.net)
        # This is just a logger warning for now
        # with pytest.raises(Exception):
        #     m.set_occupancy(phase=self.water, Pvals=self.net.Ps)
        m.set_occupancy(phase=self.water, Pvals=self.net.Ps)
        m.set_occupancy(phase=self.water, Tvals=self.net.Ts)


if __name__ == '__main__':

    t = MultiPhaseTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

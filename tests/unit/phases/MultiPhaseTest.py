import pytest
import scipy as sp
import openpnm as op
from numpy.testing import assert_allclose


class MultiPhaseTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.water = op.phases.Water(network=self.net, name="water")
        self.air = op.phases.Air(network=self.net, name="air")

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

    def test_multiphase_set_occupancy_w_indices_only(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.air,
                                                           self.water])
        Ps_water = sp.array([0, 1, 2])
        Ps_water_mask = self.net.tomask(pores=Ps_water)
        Ts_water = sp.array([4, 12, 22])
        Ts_water_mask = self.net.tomask(throats=Ts_water)
        m.set_occupancy(self.water, pores=Ps_water, throats=Ts_water)
        assert m["pore.occupancy.water"][Ps_water_mask].mean() == 1
        assert m["pore.occupancy.water"][~Ps_water_mask].mean() == 0
        assert m["throat.occupancy.water"][Ts_water_mask].mean() == 1
        assert m["throat.occupancy.water"][~Ts_water_mask].mean() == 0

    def test_multiphase_set_occupancy_w_values_only(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.air,
                                                           self.water])
        Pvals = sp.array([0.5, 0.9, 0.01])
        Tvals = sp.array([0.9, 0.01])
        # Pvals must be Np-long if not accompanied by "pores" argument
        with pytest.raises(Exception):
            m.set_occupancy(self.water, Pvals=Pvals)
        # Tvals must be Nt-long if not accompanied by "throats" argument
        with pytest.raises(Exception):
            m.set_occupancy(self.water, Tvals=Tvals)
        # Set pore occupancy
        m.set_occupancy(self.water, Pvals=1)
        assert m["pore.occupancy.water"].mean() == 1
        Pvals = sp.ones(self.net.Np) * 0.5
        m.set_occupancy(self.water, Pvals=Pvals)
        assert m["pore.occupancy.water"].mean() == 0.5
        # Setting throat occupancy
        m.set_occupancy(self.water, Tvals=1)
        assert m["throat.occupancy.water"].mean() == 1
        Tvals = sp.ones(self.net.Nt) * 0.54
        m.set_occupancy(self.water, Tvals=Tvals)
        assert m["throat.occupancy.water"].mean() == 0.54

    def test_multiphase_set_occupancy_w_pore_indices_and_Pvals(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.air,
                                                           self.water])
        Ps_water = sp.array([0, 1, 2])
        Pvals = sp.array([0.5, 0.9, 0.01])
        Ps_water_mask = self.net.tomask(Ps_water)
        Ts_water = sp.array([4, 12, 22])
        Ts_water_mask = self.net.tomask(throats=Ts_water)
        Tvals = sp.array([0.3, 0.4, 0.1])
        # Pvals/Tvals and pores/throats; same array length
        m.set_occupancy(self.water, pores=Ps_water, Pvals=Pvals,
                        throats=Ts_water, Tvals=Tvals)
        assert_allclose(m["pore.occupancy.water"][Ps_water], Pvals)
        assert_allclose(m["throat.occupancy.water"][Ts_water], Tvals)
        assert m["pore.occupancy.water"][~Ps_water_mask].mean() == 0
        assert m["throat.occupancy.water"][~Ts_water_mask].mean() == 0
        # Pvals and pores; inconsistent size
        with pytest.raises(Exception):
            m.set_occupancy(self.water, pores=[1, 5, 10], Pvals=[0.5, 0])
        # Tvals and throats; inconsistent size
        with pytest.raises(Exception):
            m.set_occupancy(self.water, throats=[10, 52, 0], Tvals=[0.5, 0.01])
        # Pvals/Tvals.size = 1 and pores/throats.size > 1
        m.set_occupancy(self.water, pores=Ps_water, Pvals=0.25,
                        throats=Ts_water, Tvals=0.23)
        assert m["pore.occupancy.water"][Ps_water].mean() == 0.25
        assert m["throat.occupancy.water"][Ts_water].mean() == 0.23

    def test_multiphase_automatic_throat_occupancy(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.air,
                                                           self.water])
        pores = sp.random.choice(self.net.Ps, size=100, replace=False)
        Pvals = sp.random.random(pores.size)
        m.set_occupancy(self.water, pores=pores, Pvals=Pvals)
        P1, P2 = self.net["throat.conns"].T
        oc1, oc2 = [m["pore.occupancy.water"][x] for x in (P1, P2)]
        # Throats take average occupancy of adjacent pores
        m.set_automatic_throat_occupancy(mode="mean")
        assert_allclose(m["throat.occupancy.water"], (oc1+oc2)/2)
        # Throats take maximum occupancy of adjacent pores
        m.set_automatic_throat_occupancy(mode="max")
        assert_allclose(m["throat.occupancy.water"], sp.maximum(oc1, oc2))
        # Throats take minimum occupancy of adjacent pores
        m.set_automatic_throat_occupancy(mode="min")
        assert_allclose(m["throat.occupancy.water"], sp.minimum(oc1, oc2))

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

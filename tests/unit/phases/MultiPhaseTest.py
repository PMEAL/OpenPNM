import pytest
import numpy as np
import scipy as sp
import openpnm as op
from numpy.testing import assert_allclose


class MultiPhaseTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.water = op.phases.Water(network=self.net, name="water")
        self.air = op.phases.Air(network=self.net, name="air")
        self.oil = op.phases.Water(network=self.net, name="oil")

    def test_multiphase_init(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.air,
                                                           self.water])
        assert np.all(m['pore.occupancy.all'] == 0.0)
        assert np.all(m['throat.occupancy.all'] == 0.0)
        assert self.air.name in m.settings['phases']
        assert self.water.name in m.settings['phases']

    def test_multiphase_no_occupancy_yet(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.air,
                                                           self.water])
        self.water['pore.temperature'] = 300
        assert np.all(self.water['pore.temperature'] == 300)
        with pytest.raises(Exception):
            m['pore.temperature']

    def test_multiphase_set_occupancy_w_indices_only(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.air,
                                                           self.water])
        Ps_water = np.array([0, 1, 2])
        Ps_water_mask = self.net.tomask(pores=Ps_water)
        Ts_water = np.array([4, 12, 22])
        Ts_water_mask = self.net.tomask(throats=Ts_water)
        m.set_occupancy(self.water, pores=Ps_water, throats=Ts_water)
        assert m["pore.occupancy.water"][Ps_water_mask].mean() == 1
        assert m["pore.occupancy.water"][~Ps_water_mask].mean() == 0
        assert m["throat.occupancy.water"][Ts_water_mask].mean() == 1
        assert m["throat.occupancy.water"][~Ts_water_mask].mean() == 0

    def test_multiphase_set_occupancy_w_values_only(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.air,
                                                           self.water])
        Pvals = np.array([0.5, 0.9, 0.01])
        Tvals = np.array([0.9, 0.01])
        # Pvals must be Np-long if not accompanied by "pores" argument
        with pytest.raises(Exception):
            m.set_occupancy(self.water, Pvals=Pvals)
        # Tvals must be Nt-long if not accompanied by "throats" argument
        with pytest.raises(Exception):
            m.set_occupancy(self.water, Tvals=Tvals)
        # Set pore occupancy
        m.set_occupancy(self.water, Pvals=1)
        assert m["pore.occupancy.water"].mean() == 1
        Pvals = np.ones(self.net.Np) * 0.5
        m.set_occupancy(self.water, Pvals=Pvals)
        assert m["pore.occupancy.water"].mean() == 0.5
        # Setting throat occupancy
        m.set_occupancy(self.water, Tvals=1)
        assert m["throat.occupancy.water"].mean() == 1
        Tvals = np.ones(self.net.Nt) * 0.54
        m.set_occupancy(self.water, Tvals=Tvals)
        assert m["throat.occupancy.water"].mean() == 0.54

    def test_multiphase_set_occupancy_w_pore_indices_and_Pvals(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.air,
                                                           self.water])
        Ps_water = np.array([0, 1, 2])
        Pvals = np.array([0.5, 0.9, 0.01])
        Ps_water_mask = self.net.tomask(Ps_water)
        Ts_water = np.array([4, 12, 22])
        Ts_water_mask = self.net.tomask(throats=Ts_water)
        Tvals = np.array([0.3, 0.4, 0.1])
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
        m._set_automatic_throat_occupancy(mode="mean")
        assert_allclose(m["throat.occupancy.water"], (oc1+oc2)/2)
        # Throats take maximum occupancy of adjacent pores
        m._set_automatic_throat_occupancy(mode="max")
        assert_allclose(m["throat.occupancy.water"], np.maximum(oc1, oc2))
        # Throats take minimum occupancy of adjacent pores
        m._set_automatic_throat_occupancy(mode="min")
        assert_allclose(m["throat.occupancy.water"], np.minimum(oc1, oc2))

    def test_multiphase_occupancy_set_single_phase(self):
        m = op.phases.MultiPhase(network=self.net)
        self.water['pore.temperature'] = 300
        assert np.all(self.water['pore.temperature'] == 300)
        m.set_occupancy(phase=self.water, Pvals=1, Tvals=1)
        assert np.all(m['pore.temperature'] == 300)

    def test_multiphase_occupancy_set_two_phase(self):
        m = op.phases.MultiPhase(network=self.net)
        self.water['pore.temperature'] = 300
        self.air['pore.temperature'] = 200
        assert np.all(self.water['pore.temperature'] == 300)
        assert np.all(self.air['pore.temperature'] == 200)
        Ps = self.net['pore.coords'][:, 0] < 3
        Ts = self.net.tomask(throats=self.net.find_neighbor_throats(Ps))
        m.set_occupancy(phase=self.water, Pvals=Ps, Tvals=Ts)
        m.set_occupancy(phase=self.air, Pvals=~Ps, Tvals=~Ts)
        assert np.all(m['pore.temperature'] >= 200)
        assert np.all(m['pore.temperature'] <= 300)

    def test_mutliphase_partition_coef(self):
        m = op.phases.MultiPhase(network=self.net,
                                 phases=[self.water, self.air, self.oil])
        x, y, z = self.net["pore.coords"].T
        ps_water = self.net.Ps[(y <= 3) + (y >= 8)]
        ps_air = self.net.Ps[(y > 3) * (y < 6)]
        ps_oil = self.net.Ps[(y >= 6) * (y < 8)]
        # Phase arrangement (y-axis): W | A | O | W
        m.set_occupancy(phase=self.water, pores=ps_water)
        m.set_occupancy(phase=self.air, pores=ps_air)
        m.set_occupancy(phase=self.oil, pores=ps_oil)
        const = op.models.misc.constant
        K_air_water = 2.0
        K_air_oil = 1.8
        K_water_oil = 0.73
        m.set_binary_partition_coef(propname="throat.partition_coef",
                                    phases=[self.air, self.water],
                                    model=const, value=K_air_water)
        m.set_binary_partition_coef(propname="throat.partition_coef",
                                    phases=[self.air, self.oil],
                                    model=const, value=K_air_oil)
        m.set_binary_partition_coef(propname="throat.partition_coef",
                                    phases=[self.water, self.oil],
                                    model=const, value=K_water_oil)
        K_aw = m["throat.partition_coef.air:water"]
        K_ao = m["throat.partition_coef.air:oil"]
        K_wo = m["throat.partition_coef.water:oil"]
        K_global = m["throat.partition_coef.all"]
        assert np.isclose(K_aw.mean(), K_air_water)
        assert np.isclose(K_ao.mean(), K_air_oil)
        assert np.isclose(K_wo.mean(), K_water_oil)
        # Get water-air interface throats
        tmp1 = self.net.find_neighbor_throats(ps_water, mode="xor")
        tmp2 = self.net.find_neighbor_throats(ps_air, mode="xor")
        Ts_water_air_interface = np.intersect1d(tmp1, tmp2)
        # Get air-oil interface throats
        tmp1 = self.net.find_neighbor_throats(ps_air, mode="xor")
        tmp2 = self.net.find_neighbor_throats(ps_oil, mode="xor")
        Ts_air_oil_interface = np.intersect1d(tmp1, tmp2)
        # Get oil-water interface throats
        tmp1 = self.net.find_neighbor_throats(ps_oil, mode="xor")
        tmp2 = self.net.find_neighbor_throats(ps_water, mode="xor")
        Ts_oil_water_interface = np.intersect1d(tmp1, tmp2)
        # K_global for water-air interface must be 1/K_air_water
        assert np.isclose(K_global[Ts_water_air_interface].mean(), 1/K_air_water)
        # K_global for air-oil interface must be K_air_oil (not 1/K_air_oil)
        assert np.isclose(K_global[Ts_air_oil_interface].mean(), K_air_oil)
        # K_global for oil-water interface must be 1/K_water_oil
        assert np.isclose(K_global[Ts_oil_water_interface].mean(), 1/K_water_oil)
        # K_global for single-phase regions must be 1.0
        interface_throats = np.hstack((Ts_water_air_interface,
                                       Ts_air_oil_interface,
                                       Ts_oil_water_interface))
        Ts_single_phase = np.setdiff1d(self.net.Ts, interface_throats)
        assert np.isclose(K_global[Ts_single_phase].mean(), 1.0)

    def test_multiphase_invalid_phase(self):
        pn = op.network.Cubic(shape=[3, 3, 3])
        water = op.phases.Water(network=pn)
        m = op.phases.MultiPhase(network=self.net)
        with pytest.raises(Exception):
            m.set_occupancy(phase=water)

    def test_multiphase_invalid_occupancy(self):
        m = op.phases.MultiPhase(network=self.net, phases=[self.water, self.air])
        # The next line ideally should throw an Exception, but warning for now
        m.set_occupancy(phase=self.water, Pvals=1.5, Tvals=2.5)


if __name__ == '__main__':

    t = MultiPhaseTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

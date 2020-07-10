import numpy as np
import openpnm as op
from numpy.testing import assert_allclose


class DispersiveConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[2, 2, 2])
        self.geo = op.geometry.GenericGeometry(network=self.net)
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase["pore.pressure"] = np.arange(self.net.Np, dtype=float)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys["throat.diffusive_conductance"] = np.arange(self.net.Nt) + 1.0
        self.phys["throat.hydraulic_conductance"] = np.arange(self.net.Nt) + 4.0

    def test_upwind(self):
        mod = op.models.physics.dispersive_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance',
                            model=mod, s_scheme="upwind")
        actual = self.phys['throat.ad_dif_conductance'].mean()
        assert_allclose(actual, desired=19.583333)

    def test_hybrid(self):
        mod = op.models.physics.dispersive_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance',
                            model=mod, s_scheme="hybrid")
        actual = self.phys['throat.ad_dif_conductance'].mean()
        assert_allclose(actual, desired=13.125)

    def test_powerlaw(self):
        mod = op.models.physics.dispersive_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance',
                            model=mod, s_scheme="powerlaw")
        actual = self.phys['throat.ad_dif_conductance'].mean()
        assert_allclose(actual, desired=13.820509)

    def test_exponential(self):
        mod = op.models.physics.dispersive_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance',
                            model=mod, s_scheme="exponential")
        actual = self.phys['throat.ad_dif_conductance'].mean()
        assert_allclose(actual, desired=13.796022)

    def test_consistency_w_diffusive_conductance_when_no_flow(self):
        self.phase["pore.uniform_pressure"] = 0.0
        mod = op.models.physics.dispersive_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance',
                            model=mod, s_scheme="powerlaw",
                            pore_pressure="pore.uniform_pressure")
        g_ad_dif = self.phys['throat.ad_dif_conductance'].mean()
        g_dif = self.phys['throat.diffusive_conductance'].mean()
        assert_allclose(actual=g_ad_dif, desired=g_dif)

    def test_consistency_w_Nt_by_2_diffusive_conductance(self):
        gd = self.phase["throat.diffusive_conductance"]
        self.phys["throat.Nt_by_2"] = np.vstack((gd, gd)).T
        mod = op.models.physics.dispersive_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance',
                            model=mod, s_scheme="powerlaw",
                            pore_pressure="pore.uniform_pressure",
                            throat_diffusive_conductance="throat.Nt_by_2")
        g_ad_dif = self.phys['throat.ad_dif_conductance']
        g_dif = self.phys['throat.Nt_by_2']
        assert_allclose(actual=g_ad_dif, desired=g_dif)

    def test_consistency_w_assymetric_Nt_by_2_diffusive_conductance(self):
        gd = self.phase["throat.diffusive_conductance"]
        self.phys["throat.Nt_by_2"] = np.vstack((gd, gd)).T
        self.phys["throat.Nt_by_2"][3, 1] = 7.45
        mod = op.models.physics.dispersive_conductance.ad_dif
        self.phys.add_model(propname='throat.ad_dif_conductance',
                            model=mod, s_scheme="powerlaw",
                            pore_pressure="pore.uniform_pressure",
                            throat_diffusive_conductance="throat.Nt_by_2")
        g_ad_dif = self.phys['throat.ad_dif_conductance']
        g_dif = self.phys['throat.Nt_by_2']
        assert_allclose(actual=g_ad_dif, desired=g_dif)


if __name__ == '__main__':

    t = DispersiveConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

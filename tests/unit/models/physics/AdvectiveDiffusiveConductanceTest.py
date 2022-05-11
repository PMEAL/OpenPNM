import pytest
import numpy as np
import openpnm as op
from numpy.testing import assert_allclose


class Ad_diff_ConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[2, 2, 2])
        self.net['throat.conduit_lengths.pore1'] = 0.15
        self.net['throat.conduit_lengths.throat'] = 0.6
        self.net['throat.conduit_lengths.pore2'] = 0.25
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase["pore.pressure"] = np.arange(self.net.Np, dtype=float)
        self.phase["throat.diffusive_conductance"] = np.arange(self.net.Nt) + 1.0
        self.phase["throat.hydraulic_conductance"] = np.arange(self.net.Nt) + 4.0

    def test_unsupported_discretization_scheme(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        with pytest.raises(Exception):
            self.phase.add_model(propname='throat.ad_dif_conductance',
                                 model=mod, s_scheme="unsupported_scheme")

    def test_upwind(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance',
                             model=mod, s_scheme="upwind")
        actual = self.phase['throat.ad_dif_conductance'].mean()
        assert_allclose(actual, desired=19.583333)

    def test_hybrid(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance',
                             model=mod, s_scheme="hybrid")
        actual = self.phase['throat.ad_dif_conductance'].mean()
        assert_allclose(actual, desired=13.125)

    def test_powerlaw(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance',
                             model=mod, s_scheme="powerlaw")
        actual = self.phase['throat.ad_dif_conductance'].mean()
        assert_allclose(actual, desired=13.820509)

    def test_exponential(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance',
                             model=mod, s_scheme="exponential")
        actual = self.phase['throat.ad_dif_conductance'].mean()
        assert_allclose(actual, desired=13.796022)

    def test_consistency_w_diffusive_conductance_when_no_flow(self):
        self.phase["pore.uniform_pressure"] = 0.0
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance',
                             model=mod, s_scheme="powerlaw",
                             pore_pressure="pore.uniform_pressure")
        g_ad_dif = self.phase['throat.ad_dif_conductance'].mean()
        g_dif = self.phase['throat.diffusive_conductance'].mean()
        assert_allclose(actual=g_ad_dif, desired=g_dif)

    def test_consistency_w_Nt_by_2_diffusive_conductance(self):
        gd = self.phase["throat.diffusive_conductance"]
        self.phase["throat.Nt_by_2"] = np.vstack((gd, gd)).T
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance',
                             model=mod, s_scheme="powerlaw",
                             pore_pressure="pore.uniform_pressure",
                             throat_diffusive_conductance="throat.Nt_by_2")
        g_ad_dif = self.phase['throat.ad_dif_conductance']
        g_dif = self.phase['throat.Nt_by_2']
        assert_allclose(actual=g_ad_dif, desired=g_dif)

    def test_consistency_w_assymetric_Nt_by_2_diffusive_conductance(self):
        gd = self.phase["throat.diffusive_conductance"]
        self.phase["throat.Nt_by_2"] = np.vstack((gd, gd)).T
        self.phase["throat.Nt_by_2"][3, 1] = 7.45
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance',
                             model=mod, s_scheme="powerlaw",
                             pore_pressure="pore.uniform_pressure",
                             throat_diffusive_conductance="throat.Nt_by_2")
        g_ad_dif = self.phase['throat.ad_dif_conductance']
        g_dif = self.phase['throat.Nt_by_2']
        assert_allclose(actual=g_ad_dif, desired=g_dif)

    def test_ad_dif_when_dif_cond_in_wrong_shape(self):
        gd = self.phase["throat.diffusive_conductance"]
        self.phase["throat.Nt_by_3"] = np.vstack((gd, gd, gd)).T
        mod = op.models.physics.ad_dif_conductance.ad_dif
        with pytest.raises(Exception):
            self.phase.add_model(propname='throat.ad_dif_conductance',
                                 model=mod, s_scheme="powerlaw",
                                 throat_diffusive_conductance="throat.Nt_by_3")


if __name__ == '__main__':

    t = Ad_diff_ConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

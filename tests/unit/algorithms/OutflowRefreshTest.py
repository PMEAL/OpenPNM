import numpy as np
import pytest
from numpy.testing import assert_allclose

import openpnm as op


class OutflowRefreshTest:

    def setup_method(self):
        self.net = op.network.Cubic(shape=[3, 2, 1])
        self.net['pore.volume'] = 1.0
        self.phase = op.phase.Phase(network=self.net)
        self.phase['throat.hydraulic_conductance'] = 1.0
        self.phase['throat.diffusive_conductance'] = 1.0
        self.flow = op.algorithms.StokesFlow(network=self.net, phase=self.phase)
        self.inlet = self.net.pores('right')
        self.outlet = self.net.pores('left')
        self.flow.set_value_BC(pores=self.inlet, values=1.0)
        self.flow.set_value_BC(pores=self.outlet, values=0.0)
        self.flow.run()
        self.phase.update(self.flow.soln)

    def teardown_method(self):
        op.Workspace().close_project(self.net.project)

    def make_transport(self, scheme='powerlaw', transient=False):
        cls = (op.algorithms.TransientAdvectionDiffusion if transient
               else op.algorithms.AdvectionDiffusion)
        ad = cls(network=self.net, phase=self.phase)
        self.phase.add_model(
            propname='throat.ad_dif_conductance',
            model=op.models.physics.ad_dif_conductance.ad_dif,
            s_scheme=scheme,
        )
        ad.set_value_BC(pores=self.inlet, values=2.0)
        ad.set_outflow_BC(pores=self.outlet)
        return ad

    @pytest.mark.parametrize('scheme', ['upwind', 'hybrid', 'powerlaw', 'exponential'])
    @pytest.mark.parametrize('change', ['pressure', 'conductance', 'zero_flow'])
    def test_outflow_refreshes_on_rerun(self, scheme, change):
        ad = self.make_transport(scheme=scheme)
        ad.run()
        assert_allclose(ad.x, 2.0)
        if change == 'pressure':
            self.flow.set_value_BC(pores=self.inlet, values=3.0, mode='overwrite')
            expected = 1.5
        elif change == 'conductance':
            self.phase['throat.hydraulic_conductance'] *= 3.0
            self.flow.settings.cache = False
            expected = 1.5
        else:
            self.flow.set_value_BC(pores=self.inlet, values=0.0, mode='overwrite')
            expected = 0.0
        self.flow.run()
        self.phase.update(self.flow.soln)
        self.phase.regenerate_models()
        ad.run()
        # Zero outlet gradient preserves a uniform, source-free concentration.
        assert_allclose(ad.x, 2.0)
        assert_allclose(ad['pore.bc.outflow'][self.outlet], expected, atol=1e-14)
        assert np.isnan(ad['pore.bc.outflow'][self.net.pores('left', mode='not')]).all()
        # A fresh instance, with the BC set after flow changes, must agree.
        fresh = self.make_transport(scheme=scheme)
        fresh.run()
        assert_allclose(ad.x, fresh.x)

    def test_outflow_refreshes_during_transient_rhs(self):
        ad = self.make_transport(transient=True)
        rhs = ad._build_rhs()
        concentration = np.full(self.net.Np, 2.0)
        assert_allclose(rhs(0.0, concentration), 0.0, atol=1e-14)
        self.phase['pore.pressure'] *= 3.0
        self.phase.regenerate_models()
        # The same RHS closure must use the updated flow on the next call.
        assert_allclose(rhs(1.0, concentration), 0.0, atol=1e-14)
        assert_allclose(ad['pore.bc.outflow'][self.outlet], 1.5)

    def test_removed_outflow_is_not_restored(self):
        ad = self.make_transport()
        ad.set_outflow_BC(pores=self.outlet, mode='remove')
        ad.set_value_BC(pores=self.outlet, values=0.0)
        self.phase['pore.pressure'] *= 3.0
        self.phase.regenerate_models()
        ad.run()
        assert np.isnan(ad['pore.bc.outflow']).all()
        assert_allclose(ad.x[self.outlet], 0.0)

    def test_zero_outflow_remains_selected_when_flow_restarts(self):
        self.phase['pore.pressure'] = 0.0
        ad = self.make_transport()
        ad.run()
        assert_allclose(ad['pore.bc.outflow'][self.outlet], 0.0)
        self.phase.update(self.flow.soln)
        self.phase.regenerate_models()
        ad.run()
        assert_allclose(ad['pore.bc.outflow'][self.outlet], 0.5)
        assert_allclose(ad.x, 2.0)

    def test_single_throat_and_custom_flow_property_names(self):
        net = op.network.Cubic(shape=[2, 1, 1])
        phase = op.phase.Phase(network=net)
        phase['pore.custom_pressure'] = [1.0, 0.0]
        phase['throat.custom_hydraulic_conductance'] = 1.0
        ad = op.algorithms.AdvectionDiffusion(network=net, phase=phase)
        ad.settings.pressure = 'pore.custom_pressure'
        ad.settings.hydraulic_conductance = 'throat.custom_hydraulic_conductance'
        ad.set_outflow_BC(pores=[1])
        assert_allclose(ad['pore.bc.outflow'][1], 1.0)
        phase['pore.custom_pressure'] = [3.0, 0.0]
        phase['throat.ad_dif_conductance'] = 1.0
        ad._build_A()
        ad._build_b()
        ad._apply_BCs()
        assert_allclose(ad['pore.bc.outflow'][1], 3.0)
        op.Workspace().close_project(net.project)

    def test_no_outflow_does_not_require_flow_properties(self):
        ad = op.algorithms.AdvectionDiffusion(network=self.net, phase=self.phase)
        self.phase['throat.ad_dif_conductance'] = 1.0
        del self.phase['pore.pressure']
        del self.phase['throat.hydraulic_conductance']
        ad.set_value_BC(pores=self.inlet, values=2.0)
        ad.set_value_BC(pores=self.outlet, values=0.0)
        ad.run()
        assert_allclose(ad.x[self.inlet], 2.0)
        assert_allclose(ad.x[self.outlet], 0.0)

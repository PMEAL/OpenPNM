import numpy as np
import openpnm as op
from numpy.testing import assert_allclose
from openpnm.algorithms import AdvectionDiffusion, StokesFlow
import pytest


class AdvectionDiffusionTest:

    def setup_class(self):
        np.random.seed(0)
        self.net = op.network.Cubic(shape=[4, 3, 1], spacing=1.)
        self.net.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders()
        )
        self.net.regenerate_models()

        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['throat.diffusive_conductance'] = 1e-15
        self.phase['throat.hydraulic_conductance'] = 1e-15

        self.sf = StokesFlow(network=self.net, phase=self.phase)
        self.sf.set_value_BC(pores=self.net.pores('right'), values=1)
        self.sf.set_value_BC(pores=self.net.pores('left'), values=0)
        self.sf.run()

        self.phase.update(self.sf.soln)

        self.ad = AdvectionDiffusion(network=self.net, phase=self.phase)
        self.ad.set_value_BC(pores=self.net.pores('right'), values=2)
        self.ad.set_value_BC(pores=self.net.pores('left'), values=0)
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance',
                             model=mod, s_scheme='powerlaw')
        self.ad.run()

    def test_settings(self):
        self.ad.settings._update({'quantity': "pore.blah",
                                  'conductance': "throat.foo",
                                  'diffusive_conductance': "throat.foo2",
                                  'hydraulic_conductance': "throat.foo3",
                                  'pressure': "pore.bar",
                                  })
        assert self.ad.settings["quantity"] == "pore.blah"
        assert self.ad.settings["conductance"] == "throat.foo"
        assert self.ad.settings["diffusive_conductance"] == "throat.foo2"
        assert self.ad.settings["hydraulic_conductance"] == "throat.foo3"
        assert self.ad.settings["pressure"] == "pore.bar"
        # Reset back to previous values for other tests to run
        self.ad.settings._update({'quantity': "pore.concentration",
                                  'conductance': "throat.ad_dif_conductance",
                                  "diffusive_conductance": "throat.diffusive_conductance",
                                  "hydraulic_conductance": "throat.hydraulic_conductance",
                                  "pressure": "pore.pressure"
                                  })

    def test_conductance_gets_updated_when_pressure_changes(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance',
                             model=mod, s_scheme='powerlaw')
        g_old = np.copy(self.phase["throat.ad_dif_conductance"])
        # Run StokesFlow with a different BC to change pressure field
        self.sf.set_value_BC(pores=self.net.pores('right'), values=1.5)
        # Running the next line should update "throat.ad_dif_conductance"
        P = self.sf.run()
        self.phase.update(P)
        self.ad.settings["conductance"] = "throat.ad_dif_conductance"
        _ = self.ad.run()
        self.phase.regenerate_models()
        g_updated = self.phase["throat.ad_dif_conductance"]
        # Ensure conductance values are updated
        assert g_old.mean() != g_updated.mean()
        assert_allclose(g_updated.mean(), 1.01258990e-15)
        # Reset BCs for other tests to run properly
        self.sf.set_value_BC(pores=self.net.pores('right'), values=1)
        self.sf.run()

    def test_powerlaw_advection_diffusion(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance_powerlaw',
                             model=mod, s_scheme='powerlaw')
        self.phase.regenerate_models()
        self.ad.settings['conductance'] = 'throat.ad_dif_conductance_powerlaw'
        self.ad.run()
        x = [
            0., 0., 0.,
            0.89653, 0.89653, 0.89653,
            1.53924, 1.53924, 1.53924,
            2., 2., 2.
        ]
        y = self.ad['pore.concentration']
        assert_allclose(actual=y, desired=x, rtol=1e-5)

    def test_upwind_advection_diffusion(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance_upwind',
                             model=mod, s_scheme='upwind')
        self.phase.regenerate_models()
        self.ad.settings['conductance'] = 'throat.ad_dif_conductance_upwind'
        self.ad.run()
        x = [
            0., 0., 0.,
            0.86486, 0.86486, 0.86486,
            1.51351, 1.51351, 1.51351,
            2., 2., 2.
        ]
        y = self.ad['pore.concentration']
        assert_allclose(actual=y, desired=x, rtol=1e-5)

    def test_hybrid_advection_diffusion(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance_hybrid',
                             model=mod, s_scheme='hybrid')
        self.phase.regenerate_models()

        self.ad.settings['conductance'] = 'throat.ad_dif_conductance_hybrid'
        self.ad.run()
        x = [
            0., 0., 0.,
            0.89908, 0.89908, 0.89908,
            1.54128, 1.54128, 1.54128,
            2., 2., 2.
        ]
        y = self.ad['pore.concentration']
        assert_allclose(actual=y, desired=x, rtol=1e-5)

    def test_exponential_advection_diffusion(self):
        mod = op.models.physics.ad_dif_conductance.ad_dif
        self.phase.add_model(propname='throat.ad_dif_conductance_exponential',
                             model=mod, s_scheme='exponential')
        self.phase.regenerate_models()

        self.ad.settings['conductance'] = 'throat.ad_dif_conductance_exponential'
        self.ad.run()
        x = [
            0., 0., 0.,
            0.89688173, 0.89688173, 0.89688173,
            1.53952557, 1.53952557, 1.53952557,
            2., 2., 2.
        ]
        y = self.ad['pore.concentration']
        assert_allclose(actual=y, desired=x, rtol=1e-5)

    def test_outflow_BC(self):
        ad = AdvectionDiffusion(network=self.net, phase=self.phase)
        ad.settings["cache"] = False
        ad.set_value_BC(pores=self.net.pores('right'), values=2)

        for s_scheme in ['upwind', 'hybrid', 'powerlaw', 'exponential']:
            ad.settings._update({'quantity': 'pore.concentration',
                                 'conductance': f'throat.ad_dif_conductance_{s_scheme}'
            })
            ad.set_outflow_BC(pores=self.net.pores('left'))
            ad.run()
            y = ad[ad.settings['quantity']].mean()
            assert_allclose(actual=y, desired=2, rtol=1e-5)

    def test_add_outflow_overwrite_rate_and_value_BC(self):
        ad = AdvectionDiffusion(network=self.net, phase=self.phase)
        ad.set_rate_BC(pores=[0, 1], total_rate=1)
        ad.set_value_BC(pores=[2, 3], values=1)
        assert np.sum(np.isfinite(ad['pore.bc_rate'])) == 2
        assert np.sum(np.isfinite(ad['pore.bc_value'])) == 2
        ad.set_outflow_BC(pores=[0, 1, 2, 3])
        assert np.sum(np.isfinite(ad['pore.bc_rate'])) == 0
        assert np.sum(np.isfinite(ad['pore.bc_value'])) == 0

    def test_value_BC_does_not_overwrite_outflow(self):
        ad = AdvectionDiffusion(network=self.net, phase=self.phase)
        ad.set_outflow_BC(pores=[0, 1])
        with pytest.raises(Exception):
            ad.set_value_BC(pores=[0, 1], values=1)

    def test_add_rate_BC_fails_when_outflow_BC_present(self):
        ad = AdvectionDiffusion(network=self.net, phase=self.phase)
        ad.set_outflow_BC(pores=[0, 1])
        with pytest.raises(Exception):
            ad.set_rate_BC(pores=[0, 1], total_rate=1)
        ad.set_rate_BC(pores=[2, 3], total_rate=1)
        assert np.all(ad['pore.bc_rate'][[2, 3]] == 0.5)

    def test_outflow_BC_rigorous(self):
        ad = AdvectionDiffusion(network=self.net, phase=self.phase)
        ad.settings["cache"] = False
        # Add source term so we get a non-uniform concentration profile
        self.phase["pore.A1"] = -5e-16
        linear = op.models.physics.generic_source_term.linear
        self.phase.add_model(
            propname="pore.rxn",
            model=linear,
            A1="pore.A1",
            X="pore.concentration",
            regen_mode='deferred',
        )
        internal_pores = self.net.pores(["left", "right"], mode="not")
        ad.set_source("pore.rxn", pores=internal_pores)
        ad.set_value_BC(pores=self.net.pores('right'), values=2)
        ad.set_outflow_BC(pores=self.net.pores('left'))

        for s_scheme in ['upwind', 'hybrid', 'powerlaw', 'exponential']:
            ad.settings._update({'quantity': 'pore.concentration',
                                 'conductance': f'throat.ad_dif_conductance_{s_scheme}'})
            ad.run()
            y = ad[ad.settings['quantity']].reshape(4, 3).mean(axis=1)
            # Calculate grad_c @ face on which outflow BC was imposed
            dydx_left_mean = (y[1] - y[0]) / (1 - 0)
            # Calculate grad_c across the entire network as reference
            dydx_total_mean = (y[1] - y[-1]) / (3 - 0)
            # Make sure that grad_c @ outflow BC face is "numerically" ~ 0
            assert_allclose(
                actual=dydx_total_mean + dydx_left_mean,
                desired=dydx_total_mean
            )

    def test_rate(self):
        ad = AdvectionDiffusion(network=self.net, phase=self.phase)
        ad.settings["cache"] = False
        ad.set_value_BC(pores=self.net.pores('right'), values=2)
        ad.set_value_BC(pores=self.net.pores('left'), values=0)

        for s_scheme in ['upwind', 'hybrid', 'powerlaw', 'exponential']:
            ad.settings._update({'quantity': 'pore.concentration',
                                 'conductance': f'throat.ad_dif_conductance_{s_scheme}',
                                 's_scheme': s_scheme})
            ad.run()

            mdot_inlet = ad.rate(pores=self.net.pores("right"))[0]
            mdot_outlet = ad.rate(pores=self.net.pores("left"))[0]
            temp = np.random.choice(
                self.net.pores(["right", "left"], mode="not"), size=3, replace=False
            )
            mdot_internal = ad.rate(pores=temp)[0]
            # Ensure no mass is generated within the network
            assert_allclose(mdot_inlet - mdot_internal, mdot_inlet)
            # Ensure mass is not conserved
            assert_allclose(mdot_inlet, -mdot_outlet)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':
    t = AdvectionDiffusionTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()

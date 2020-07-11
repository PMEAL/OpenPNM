import numpy as np
import numpy.testing as nt
import openpnm as op
from openpnm.phases import mixtures


class TransientIonicTransportTest:
    def setup_class(self):
        np.random.seed(0)
        self.net = op.network.Cubic(shape=[8, 8, 1], spacing=9e-4)

        prs = (
            self.net["pore.back"] * self.net["pore.right"]
            + self.net["pore.back"] * self.net["pore.left"]
            + self.net["pore.front"] * self.net["pore.right"]
            + self.net["pore.front"] * self.net["pore.left"]
        )
        thrts = self.net["throat.surface"]
        op.topotools.trim(
            network=self.net, pores=self.net.Ps[prs], throats=self.net.Ts[thrts]
        )

        self.geo = op.geometry.StickAndBall(
            network=self.net, pores=self.net.Ps, throats=self.net.Ts
        )

        pr_d = op.models.misc.constant
        trt_d = op.models.misc.constant
        self.geo.add_model(propname="pore.diameter", model=pr_d, value=1.5e-4)
        self.geo.add_model(propname="throat.diameter", model=trt_d, value=1e-4)
        self.geo.regenerate_models()

        self.phase = mixtures.SalineWater(network=self.net)
        # Retrieve handles to each species for use below
        self.Na = self.phase.components['Na_' + self.phase.name]
        self.Cl = self.phase.components['Cl_' + self.phase.name]
        self.H2O = self.phase.components['H2O_' + self.phase.name]

        # physics
        self.phys = op.physics.GenericPhysics(
            network=self.net, phase=self.phase, geometry=self.geo
        )

        flow = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.phys.add_model(
            propname="throat.hydraulic_conductance",
            pore_viscosity="pore.viscosity",
            throat_viscosity="throat.viscosity",
            model=flow,
            regen_mode="normal",
        )

        current = op.models.physics.ionic_conductance.electroneutrality
        self.phys.add_model(
            propname="throat.ionic_conductance",
            ions=[self.Na.name, self.Cl.name],
            model=current,
            regen_mode="normal",
        )

        eA_dif = op.models.physics.diffusive_conductance.ordinary_diffusion
        self.phys.add_model(
            propname="throat.diffusive_conductance." + self.Na.name,
            pore_diffusivity="pore.diffusivity." + self.Na.name,
            throat_diffusivity="throat.diffusivity." + self.Na.name,
            model=eA_dif,
            regen_mode="normal",
        )

        eB_dif = op.models.physics.diffusive_conductance.ordinary_diffusion
        self.phys.add_model(
            propname="throat.diffusive_conductance." + self.Cl.name,
            pore_diffusivity="pore.diffusivity." + self.Cl.name,
            throat_diffusivity="throat.diffusivity." + self.Cl.name,
            model=eB_dif,
            regen_mode="normal",
        )

        # settings for algorithms
        self.settings1 = {
            "solver_max_iter": 5,
            "solver_tol": 1e-08,
            "solver_rtol": 1e-08,
            "nlin_max_iter": 10,
        }
        self.settings2 = {
            "g_tol": 1e-4,
            "g_max_iter": 4,
            "t_output": 500,
            "t_step": 500,
            "t_final": 2000,
            "t_scheme": "implicit",
        }

    def test_transient_ionic_reactive_transport(self):
        # algorithms
        sf = op.algorithms.StokesFlow(
            network=self.net, phase=self.phase, settings=self.settings1
        )
        sf.set_value_BC(pores=self.net.pores("back"), values=0.01)
        sf.set_value_BC(pores=self.net.pores("front"), values=0.00)
        sf.run()
        self.phase.update(sf.results())

        p = op.algorithms.TransientIonicConduction(
            network=self.net, phase=self.phase, settings=self.settings1
        )
        p.set_value_BC(pores=self.net.pores("left"), values=0.1)
        p.set_value_BC(pores=self.net.pores("right"), values=0.00)
        p.settings["charge_conservation"] = "electroneutrality"

        eA = op.algorithms.TransientNernstPlanck(
            network=self.net,
            phase=self.phase,
            ion=self.Na.name,
            settings=self.settings1,
        )
        eA.set_value_BC(pores=self.net.pores("back"), values=100)
        eA.set_value_BC(pores=self.net.pores("front"), values=90)

        eB = op.algorithms.TransientNernstPlanck(
            network=self.net,
            phase=self.phase,
            ion=self.Cl.name,
            settings=self.settings1,
        )
        eB.set_value_BC(pores=self.net.pores("back"), values=100)
        eB.set_value_BC(pores=self.net.pores("front"), values=90)

        ad_dif_mig_Na = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(
            propname="throat.ad_dif_mig_conductance." + self.Na.name,
            pore_pressure=sf.settings["quantity"],
            model=ad_dif_mig_Na,
            ion=self.Na.name,
            s_scheme="powerlaw",
        )

        ad_dif_mig_Cl = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(
            propname="throat.ad_dif_mig_conductance." + self.Cl.name,
            pore_pressure=sf.settings["quantity"],
            model=ad_dif_mig_Cl,
            ion=self.Cl.name,
            s_scheme="powerlaw",
        )

        it = op.algorithms.TransientNernstPlanckMultiphysics(
            network=self.net, phase=self.phase, settings=self.settings2
        )
        it.setup(potential_field=p.name, ions=[eA.name, eB.name])
        it.run()

        self.phase.update(sf.results())
        self.phase.update(p.results())
        self.phase.update(eA.results())
        self.phase.update(eB.results())

        y1 = p[p.settings["quantity"]].mean()
        y2 = eA[eA.settings["quantity"]].mean()
        y3 = eB[eB.settings["quantity"]].mean()

        nt.assert_allclose(y1, 0.04590198)
        nt.assert_allclose(y2, 77.1751280)
        nt.assert_allclose(y3, 84.0102562)

        times = [
            "pore.concentration.Na_mix_01@500",
            "pore.concentration.Cl_mix_01@500",
            "pore.potential@500",
            "pore.concentration.Na_mix_01@1000",
            "pore.concentration.Cl_mix_01@1000",
            "pore.potential@1000",
            "pore.concentration.Na_mix_01@1500",
            "pore.concentration.Cl_mix_01@1500",
            "pore.potential@1500",
            "pore.concentration.Na_mix_01@2000",
            "pore.concentration.Cl_mix_01@2000",
            "pore.potential@2000",
        ]
        assert set(times).issubset(set(self.phase.keys()))

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == "__main__":

    t = TransientIonicTransportTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith("test"):
            print("running test: " + item)
            t.__getattribute__(item)()
    self = t

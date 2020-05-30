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
            "solver_maxiter": 5,
            "solver_tol": 1e-08,
            "solver_rtol": 1e-08,
            "max_iter": 10,
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

        x = [8.52358641e-02, 7.09614385e-02, 5.69582435e-02, 4.30417565e-02,
             2.90385615e-02, 1.47641359e-02, 1.00000000e-01, 8.52358641e-02,
             7.09614385e-02, 5.69582435e-02, 4.30417565e-02, 2.90385615e-02,
             1.47641359e-02, 0.00000000e+00, 1.00000000e-01, 8.51602277e-02,
             7.08999801e-02, 5.69362166e-02, 4.30637834e-02, 2.91000199e-02,
             1.48397723e-02, 0.00000000e+00, 1.00000000e-01, 8.51289913e-02,
             7.08745683e-02, 5.69272002e-02, 4.30727998e-02, 2.91254317e-02,
             1.48710087e-02, 0.00000000e+00, 1.00000000e-01, 8.51678142e-02,
             7.09107984e-02, 5.69412249e-02, 4.30587751e-02, 2.90892016e-02,
             1.48321858e-02, 0.00000000e+00, 1.00000000e-01, 8.52557265e-02,
             7.09902337e-02, 5.69713434e-02, 4.30286566e-02, 2.90097663e-02,
             1.47442735e-02, 0.00000000e+00, 1.00000000e-01, 8.53454249e-02,
             7.10679654e-02, 5.70001588e-02, 4.29998412e-02, 2.89320346e-02,
             1.46545751e-02, 0.00000000e+00, 8.53454249e-02, 7.10679654e-02,
             5.70001588e-02, 4.29998412e-02, 2.89320346e-02, 1.46545751e-02,
             9.00000000e+01, 9.00000000e+01, 9.00000000e+01, 9.00000000e+01,
             9.00000000e+01, 9.00000000e+01, 6.83436209e+01, 7.21420259e+01,
             7.42955208e+01, 7.52813639e+01, 7.52813639e+01, 7.42955209e+01,
             7.21420259e+01, 6.83436209e+01, 5.91341998e+01, 6.39862826e+01,
             6.70219572e+01, 6.84684599e+01, 6.84684599e+01, 6.70219572e+01,
             6.39862826e+01, 5.91341998e+01, 5.84919166e+01, 6.36762863e+01,
             6.69726467e+01, 6.85549645e+01, 6.85549645e+01, 6.69726467e+01,
             6.36762863e+01, 5.84919166e+01, 6.39543314e+01, 6.90739442e+01,
             7.22546247e+01, 7.37556768e+01, 7.37556768e+01, 7.22546247e+01,
             6.90739442e+01, 6.39543314e+01, 7.36077215e+01, 7.82086043e+01,
             8.08818628e+01, 8.20902100e+01, 8.20902100e+01, 8.08818628e+01,
             7.82086043e+01, 7.36077215e+01, 8.58271178e+01, 8.91971305e+01,
             9.08393149e+01, 9.15219751e+01, 9.15219751e+01, 9.08393149e+01,
             8.91971305e+01, 8.58271178e+01, 1.00000000e+02, 1.00000000e+02,
             1.00000000e+02, 1.00000000e+02, 1.00000000e+02, 1.00000000e+02,
             9.00000000e+01, 9.00000000e+01, 9.00000000e+01, 9.00000000e+01,
             9.00000000e+01, 9.00000000e+01, 7.83302919e+01, 8.04800087e+01,
             8.17344181e+01, 8.23177287e+01, 8.23177287e+01, 8.17344181e+01,
             8.04800087e+01, 7.83302919e+01, 7.24390267e+01, 7.54511112e+01,
             7.73402378e+01, 7.82449570e+01, 7.82449571e+01, 7.73402378e+01,
             7.54511112e+01, 7.24390267e+01, 7.16218710e+01, 7.49238264e+01,
             7.70257739e+01, 7.80401361e+01, 7.80401361e+01, 7.70257740e+01,
             7.49238264e+01, 7.16218710e+01, 7.50386393e+01, 7.82298025e+01,
             8.02364363e+01, 8.11963705e+01, 8.11963705e+01, 8.02364363e+01,
             7.82298025e+01, 7.50386393e+01, 8.16502443e+01, 8.43602521e+01,
             8.59917005e+01, 8.67513862e+01, 8.67513862e+01, 8.59917005e+01,
             8.43602521e+01, 8.16502443e+01, 9.02883066e+01, 9.20909863e+01,
             9.30507321e+01, 9.34734416e+01, 9.34734416e+01, 9.30507321e+01,
             9.20909863e+01, 9.02883066e+01, 1.00000000e+02, 1.00000000e+02,
             1.00000000e+02, 1.00000000e+02, 1.00000000e+02, 1.00000000e+02]

        y1 = p[p.settings["quantity"]]
        y2 = eA[eA.settings["quantity"]]
        y3 = eB[eB.settings["quantity"]]
        y = np.concatenate((y1, y2, y3))
        nt.assert_allclose(y, x, rtol=1e-5)

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

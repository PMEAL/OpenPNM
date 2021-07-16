import numpy as np
import openpnm as op
from openpnm.phases import mixtures
from numpy.testing import assert_allclose


class TransientMultiphysicsNernstPlanckSolverTest:

    def setup_class(self):
        # network
        np.random.seed(0)
        self.net = op.network.Cubic(shape=[6, 6, 1], spacing=1e-6)
        prs = (
            self.net['pore.back'] * self.net['pore.right']
            + self.net['pore.back'] * self.net['pore.left']
            + self.net['pore.front'] * self.net['pore.right']
            + self.net['pore.front'] * self.net['pore.left']
        )
        prs = self.net.Ps[prs]
        thrts = self.net['throat.surface']
        thrts = self.net.Ts[thrts]
        op.topotools.trim(network=self.net, pores=prs, throats=thrts)
        np.random.seed(0)
        op.topotools.reduce_coordination(self.net, 3)

        # geometry
        np.random.seed(0)
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)

        # phase
        self.sw = mixtures.SalineWater(network=self.net)
        # Retrieve handles to each species for use below
        self.Na = self.sw.components['Na_' + self.sw.name]
        self.Cl = self.sw.components['Cl_' + self.sw.name]
        self.H2O = self.sw.components['H2O_' + self.sw.name]

        # physics
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.sw,
                                              geometry=self.geo)
        modphys = op.models.physics
        self.flow = modphys.hydraulic_conductance.hagen_poiseuille
        self.phys.add_model(propname='throat.hydraulic_conductance',
                            pore_viscosity='pore.viscosity',
                            throat_viscosity='throat.viscosity',
                            model=self.flow, regen_mode='normal')
        self.current = modphys.ionic_conductance.electroneutrality
        self.phys.add_model(propname='throat.ionic_conductance',
                            ions=[self.Na.name, self.Cl.name],
                            model=self.current, regen_mode='normal')
        self.eA_dif = modphys.diffusive_conductance.ordinary_diffusion
        self.phys.add_model(propname=('throat.diffusive_conductance.' + self.Na.name),
                            pore_diffusivity=('pore.diffusivity.' + self.Na.name),
                            throat_diffusivity=('throat.diffusivity.' + self.Na.name),
                            model=self.eA_dif, regen_mode='normal')
        self.eB_dif = modphys.diffusive_conductance.ordinary_diffusion
        self.phys.add_model(propname=('throat.diffusive_conductance.' + self.Cl.name),
                            pore_diffusivity=('pore.diffusivity.' + self.Cl.name),
                            throat_diffusivity=('throat.diffusivity.' + self.Cl.name),
                            model=self.eB_dif, regen_mode='normal')

        scheme = 'powerlaw'
        self.ad_dif_mig_Na = modphys.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(propname=('throat.ad_dif_mig_conductance.' + self.Na.name),
                            pore_pressure='pore.pressure',
                            model=self.ad_dif_mig_Na, ion=self.Na.name,
                            s_scheme=scheme)

        self.ad_dif_mig_Cl = modphys.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(propname=('throat.ad_dif_mig_conductance.' + self.Cl.name),
                            pore_pressure='pore.pressure',
                            model=self.ad_dif_mig_Cl, ion=self.Cl.name,
                            s_scheme=scheme)

        # settings for algorithms
        setts1 = {'solver_max_iter': 5, 'solver_tol': 1e-08,
                  'solver_rtol': 1e-08, 'nlin_max_iter': 10,
                  'cache_A': False, 'cache_b': False}
        setts2 = {'g_tol': 1e-4, 'g_max_iter': 100, 't_output': 1000,
                  't_step': 500, 't_final': 20000, 't_scheme': 'implicit'}

        # algorithms
        self.sf = op.algorithms.StokesFlow(network=self.net, phase=self.sw,
                                           settings=setts1)
        self.sf.set_value_BC(pores=self.net.pores('right'), values=11)
        self.sf.set_value_BC(pores=self.net.pores('left'), values=10)

        self.p = op.algorithms.TransientIonicConduction(network=self.net,
                                                        phase=self.sw,
                                                        settings=setts1)
        self.p.set_value_BC(pores=self.net.pores('front'), values=0.02)
        self.p.set_value_BC(pores=self.net.pores('back'), values=0.01)
        self.p.settings['charge_conservation'] = 'laplace'

        self.eA = op.algorithms.TransientNernstPlanck(network=self.net,
                                                      phase=self.sw,
                                                      ion=self.Na.name,
                                                      settings=setts1)
        self.eA.set_value_BC(pores=self.net.pores('right'), values=20)
        self.eA.set_value_BC(pores=self.net.pores('left'), values=10)

        self.eB = op.algorithms.TransientNernstPlanck(network=self.net,
                                                      phase=self.sw,
                                                      ion=self.Cl.name,
                                                      settings=setts1)
        self.eB.set_value_BC(pores=self.net.pores('right'), values=20)
        self.eB.set_value_BC(pores=self.net.pores('left'), values=10)

        mnp = op.algorithms.TransientNernstPlanckMultiphysicsSolver
        self.mnp = mnp(network=self.net, phase=self.sw, settings=setts2)
        self.mnp.setup(potential_field=self.p.name,
                       ions=[self.eA.name, self.eB.name])

    def test_run_algs(self):
        self.sf.run()
        self.sw.update(self.sf.results())
        self.mnp.run()
        self.sw.update(self.p.results())
        self.sw.update(self.eA.results())
        self.sw.update(self.eB.results())

    def test_concentration_Na(self):
        x = [10.,         10.,         10.,         10.,        10.60812616,
             11.06857002, 11.67281073, 11.94601288, 11.8382579,  12.82984479,
             11.65799543, 12.60587744, 13.18301294, 13.44969262, 14.10341601,
             15.29681638, 13.15460857, 14.48820857, 14.48971804, 15.02531941,
             16.17694508, 17.44018265, 15.03964099, 16.23394014, 16.69406773,
             16.57725792, 17.99662053, 19.31869847, 20.,         20.,
             20.,         20.]
        x = np.around(x, decimals=5)
        y = np.around(self.sw['pore.concentration.Na_mix_01'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_concentration_Cl(self):
        x = [10.,         10.,         10.,         10.,         12.97729747,
             12.43745204, 11.61660455, 11.42491776, 10.68807667,  9.86202172,
             15.67930533, 14.50032104, 13.54875916, 12.82181397, 12.02977599,
             11.09125788, 17.79842716, 16.16013059, 14.85146661, 14.02103159,
             13.3301992,  12.36465837, 19.2899343,  17.87081166, 17.05055164,
             15.82065577, 15.89997926, 14.81186187, 20.,         20.,
             20.,         20.]
        x = np.around(x, decimals=5)
        y = np.around(self.sw['pore.concentration.Cl_mix_01'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_potential(self):
        x = [0.01890885, 0.01437119, 0.01331824, 0.01206583, 0.02,
             0.01890885, 0.01437119, 0.01331824, 0.01206583, 0.01,
             0.02,       0.01799239, 0.01545302, 0.01383603, 0.01208613,
             0.01,       0.02,       0.01751992, 0.01564312, 0.01374991,
             0.01193103, 0.01,       0.02,       0.01803751, 0.01590206,
             0.01378346, 0.01182057, 0.01,       0.01803751, 0.01590206,
             0.01378346, 0.01182057]
        x = np.around(x, decimals=5)
        y = np.around(self.sw['pore.potential'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_times(self):
        times = [
            "pore.concentration.Na_mix_01@1000",
            "pore.concentration.Cl_mix_01@1000",
            "pore.potential@1000",
            "pore.concentration.Na_mix_01@2000",
            "pore.concentration.Cl_mix_01@2000",
            "pore.potential@2000",
            "pore.concentration.Na_mix_01@3000",
            "pore.concentration.Cl_mix_01@3000",
            "pore.potential@3000",
            "pore.concentration.Na_mix_01@3500",
            "pore.concentration.Cl_mix_01@3500",
            "pore.potential@3500",
        ]
        assert set(times).issubset(set(self.sw.keys()))

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':
    t = TransientMultiphysicsNernstPlanckSolverTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

import numpy as np
import openpnm as op
from openpnm.phases import mixtures
from numpy.testing import assert_allclose


class MultiphysicsNernstPlanckSolverTest:

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
        self.phys.add_model(propname=('throat.ad_dif_mig_conductance.'
                                      + self.Cl.name),
                            pore_pressure='pore.pressure',
                            model=self.ad_dif_mig_Cl, ion=self.Cl.name,
                            s_scheme=scheme)

        # settings for algorithms
        setts1 = {'solver_max_iter': 5, 'solver_tol': 1e-08,
                  'solver_rtol': 1e-08, 'nlin_max_iter': 10,
                  'cache_A': False, 'cache_b': False}
        setts2 = {'g_tol': 1e-4, 'g_max_iter': 100}

        # algorithms
        self.sf = op.algorithms.StokesFlow(network=self.net, phase=self.sw,
                                           settings=setts1)
        self.sf.set_value_BC(pores=self.net.pores('back'), values=11)
        self.sf.set_value_BC(pores=self.net.pores('front'), values=10)

        self.p = op.algorithms.IonicConduction(network=self.net, phase=self.sw,
                                               settings=setts1)
        self.p.set_value_BC(pores=self.net.pores('left'), values=0.02)
        self.p.set_value_BC(pores=self.net.pores('right'), values=0.01)
        self.p.settings['charge_conservation'] = 'electroneutrality'

        self.eA = op.algorithms.NernstPlanck(network=self.net, phase=self.sw,
                                             ion=self.Na.name, settings=setts1)
        self.eA.set_value_BC(pores=self.net.pores('back'), values=20)
        self.eA.set_value_BC(pores=self.net.pores('front'), values=10)

        self.eB = op.algorithms.NernstPlanck(network=self.net, phase=self.sw,
                                             ion=self.Cl.name, settings=setts1)
        self.eB.set_value_BC(pores=self.net.pores('back'), values=20)
        self.eB.set_value_BC(pores=self.net.pores('front'), values=10)

        mnp = op.algorithms.NernstPlanckMultiphysicsSolver
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
        x = [16.75211024, 12.05078194, 11.35129692, 10.82954036, 20.,
             17.88816818, 12.87093078, 12.13946951, 11.38196949, 10.,
             20.,         17.48986699, 14.6700871,  13.11729069, 11.69505888,
             10.,         20.,         17.22150424, 15.25067879, 13.45710103,
             11.86509361, 10.,         20.,         18.72761149, 16.68478905,
             14.26063817, 12.41025969, 10.,         20.97223212, 18.7905309,
             16.52426576, 14.01925674]
        x = np.around(x, decimals=5)
        y = np.around(self.sw['pore.concentration.Na_mix_01'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_concentration_Cl(self):
        x = [20.83282456, 16.07633318, 14.93737703, 13.77151411, 20.,
             19.50975474, 15.05193283, 13.96754625, 13.10310732, 10.,
             20.,         17.80316316, 15.26570525, 13.86995902, 12.31635857,
             10.,         20.,         16.34464455, 14.825669,   13.18827827,
             11.65941891, 10.,         20.,         16.55949282, 14.15386805,
             12.49365792, 10.91545082, 10.,         14.78715982, 12.56772913,
             10.7821756,  9.66267911]
        x = np.around(x, decimals=5)
        y = np.around(self.sw['pore.concentration.Cl_mix_01'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_potential(self):
        x = [0.01870404, 0.01418691, 0.01324578, 0.01238089, 0.02,
             0.01870404, 0.01418691, 0.01324578, 0.01238089, 0.01,
             0.02,       0.01774593, 0.01519198, 0.0137006,  0.01212227,
             0.01,       0.02,       0.01697308, 0.01525899, 0.01349426,
             0.01185306, 0.01,       0.02,       0.01777765, 0.01560814,
             0.0135155,  0.01169786, 0.01,       0.01777765, 0.01560814,
             0.0135155,  0.01169786]
        x = np.around(x, decimals=5)
        y = np.around(self.sw['pore.potential'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':
    t = MultiphysicsNernstPlanckSolverTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

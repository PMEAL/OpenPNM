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
        self.sf.set_value_BC(pores=self.net.pores('right'), values=11)
        self.sf.set_value_BC(pores=self.net.pores('left'), values=10)

        self.p = op.algorithms.IonicConduction(network=self.net, phase=self.sw,
                                               settings=setts1)
        self.p.set_value_BC(pores=self.net.pores('front'), values=0.02)
        self.p.set_value_BC(pores=self.net.pores('back'), values=0.01)
        self.p.settings['charge_conservation'] = 'laplace'

        self.eA = op.algorithms.NernstPlanck(network=self.net, phase=self.sw,
                                             ion=self.Na.name, settings=setts1)
        self.eA.set_value_BC(pores=self.net.pores('right'), values=20)
        self.eA.set_value_BC(pores=self.net.pores('left'), values=10)

        self.eB = op.algorithms.NernstPlanck(network=self.net, phase=self.sw,
                                             ion=self.Cl.name, settings=setts1)
        self.eB.set_value_BC(pores=self.net.pores('right'), values=20)
        self.eB.set_value_BC(pores=self.net.pores('left'), values=10)

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
        x = [10.,         10.,         10.,         10.,         10.53517114,
             11.08046064, 11.71018632, 11.96714113, 11.72777708, 12.86695077,
             11.58746158, 12.65040345, 13.25574649, 13.47731388, 14.06090075,
             15.27217686, 13.05944438, 14.69280374, 14.62286844, 15.10186986,
             16.15146162, 17.35993123, 14.90573687, 16.25298948, 16.74426472,
             16.63951847, 17.98769641, 19.21709326, 20.,         20.,
             20.,         20.]
        x = np.around(x, decimals=5)
        y = np.around(self.sw['pore.concentration.Na_mix_01'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_concentration_Cl(self):
        x = [10.,         10.,         10.,         10.,         13.06578514,
             12.42279423, 11.58371717, 11.40980409, 10.79147327,  9.83605168,
             15.77521525, 14.44971313, 13.47980971, 12.80209187, 12.07204557,
             11.11458021, 17.92765688, 15.93468763, 14.72209168, 13.95745968,
             13.35854227, 12.42861968, 19.45846835, 17.84550525, 17.00086559,
             15.76790954, 15.91290826, 14.89489377, 20.,         20.,
             20.,         20.]
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
    
    def test_charge_conservation_electroneutrality(self):
        model = op.models.physics.generic_source_term.charge_conservation
        self.phys.add_model(propname='pore.charge_conservation',
                       model=model,
                       phase=self.sw,
                       p_alg=self.p,
                       e_alg=[self.eA, self.eB],
                       assumption='electroneutrality')
        x=[-2.24609335e-11, -4.86675995e-12, -2.00155305e-12,  5.17199639e-12,
          9.83787767e-12,  9.51427697e-14,  9.68585899e-14,  2.90865203e-14,
         -5.08595201e-14, -1.75206280e-11,  2.72090557e-11, -1.64924728e-13,
         -8.86829423e-15, -1.08748665e-13, -5.41010482e-14, -4.43723556e-12,
          8.16466658e-12, -4.01358941e-14, -1.27058643e-14, -5.93906159e-14,
          6.65433344e-14, -4.29943549e-11,  2.86333670e-11, -1.29223441e-13,
          2.90086871e-13, -3.21177535e-14, -1.23355336e-13, -8.68490251e-12,
         -3.95913957e-12,  3.01316467e-12,  1.28504967e-11,  1.22515955e-11]
        x = np.around(x, decimals=5)
        y = np.around(self.phys['pore.charge_conservation.rate'], decimals=5)
        assert_allclose(actual=y, desired=x)
        
    def test_charge_conservation_poisson(self):
        model = op.models.physics.generic_source_term.charge_conservation
        self.phys.add_model(propname='pore.charge_conservation',
                       model=model,
                       phase=self.sw,
                       p_alg=self.p,
                       e_alg=[self.eA, self.eB],
                       assumption='poisson')
        x=[-1.26217745e-29, -1.26217745e-29, -1.26217745e-29, -1.26217745e-29,
         -8.92959041e-15, -9.69825374e-15,  4.69293590e-16,  7.58658706e-15,
          1.49936602e-14,  9.20367258e-15, -4.47588842e-14, -9.10700568e-15,
         -1.28359063e-15,  9.93231999e-15,  1.31260953e-15,  3.03490459e-15,
         -2.28120674e-15, -1.46871791e-14, -1.02469014e-15,  1.48038702e-14,
          4.62126084e-14,  5.36985074e-14, -1.83813160e-14, -1.65456504e-14,
         -2.25583527e-16,  6.19002886e-15,  2.10179612e-15,  6.63622097e-14,
         -1.26217745e-29, -1.26217745e-29, -1.26217745e-29, -5.04870979e-29]
        x = np.around(x, decimals=5)
        y = np.around(self.phys['pore.charge_conservation.rate'], decimals=5)
        assert_allclose(actual=y, desired=x)
        self.phys.add_model(propname='pore.charge_conservation',
                       model=model,
                       phase=self.sw,
                       p_alg=self.p,
                       e_alg=[self.eA, self.eB],
                       assumption='poisson_2D')
        x=[-2.64697796e-23, -5.29395592e-23, -2.64697796e-23, -2.64697796e-23,
         -3.25242702e-08, -2.78180753e-08,  1.68087663e-09,  1.76190098e-08,
          3.29854138e-08,  3.52432404e-08, -1.12674189e-07, -2.94123872e-08,
         -3.97788700e-09,  2.24781499e-08,  8.35993025e-09,  1.86905484e-08,
         -1.62858653e-08, -3.57462500e-08, -2.60922022e-09,  3.49694375e-08,
          1.00562876e-07,  1.34340467e-07, -6.40107592e-08, -4.20463571e-08,
         -1.30577653e-09,  1.78572086e-08,  1.16045831e-08,  1.48056244e-07,
         -5.29395592e-23,  0.00000000e+00, -2.64697796e-23, -1.05879118e-22,]
        x = np.around(x, decimals=5)
        y = np.around(self.phys['pore.charge_conservation.rate'], decimals=5)
        assert_allclose(actual=y, desired=x)

if __name__ == '__main__':
    t = MultiphysicsNernstPlanckSolverTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

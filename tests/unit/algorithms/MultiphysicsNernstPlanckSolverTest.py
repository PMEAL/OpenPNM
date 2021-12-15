import numpy as np
import openpnm as op
from openpnm.phase import mixtures
from numpy.testing import assert_allclose
import pytest


@pytest.mark.skip(reason="Needs to be refactored using Integrators")
class MultiphysicsNernstPlanckSolverTest:

    def setup_class(self):
        # Create network
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

        # Create geometry
        np.random.seed(0)
        self.geo = op.geometry.SpheresAndCylinders(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)

        # Create phase
        self.sw = mixtures.SalineWater(network=self.net)
        # Retrieve handles to each species for use below
        self.Na = self.sw.components['Na_' + self.sw.name]
        self.Cl = self.sw.components['Cl_' + self.sw.name]
        self.H2O = self.sw.components['H2O_' + self.sw.name]

        # Create physics
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

        # Define settings for algorithms
        setts1 = {'solver_max_iter': 5, 'solver_tol': 1e-08,
                  'solver_rtol': 1e-08, 'nlin_max_iter': 10,
                  'cache': False}
        setts2 = {'g_tol': 1e-4, 'g_max_iter': 100}

        # Set up algorithms
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
        self.mnp.settings._update({'potential_field': self.p.name,
                                   'ions': [self.eA.name, self.eB.name]})

    def test_run_algs(self):
        self.sf.run()
        self.sw.update(self.sf.results())
        self.mnp.run()
        self.sw.update(self.p.results())
        self.sw.update(self.eA.results())
        self.sw.update(self.eB.results())

    def test_concentration_Na(self):
        y = self.sw['pore.concentration.Na_mix_01'].mean()
        assert_allclose(actual=y, desired=14.492374, rtol=1e-5)

    def test_concentration_Cl(self):
        y = self.sw['pore.concentration.Cl_mix_01'].mean()
        assert_allclose(actual=y, desired=14.335148, rtol=1e-5)

    def test_potential(self):
        y = self.sw['pore.potential'].mean()
        assert_allclose(actual=y, desired=0.0147257, rtol=1e-5)

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
        y = np.linalg.norm(self.phys['pore.charge_conservation.rate'])
        assert_allclose(actual=y, desired=6.357879e-11, rtol=1e-5)

    def test_charge_conservation_poisson(self):
        model = op.models.physics.generic_source_term.charge_conservation
        self.phys.add_model(propname='pore.charge_conservation',
                            model=model,
                            phase=self.sw,
                            p_alg=self.p,
                            e_alg=[self.eA, self.eB],
                            assumption='poisson')
        y = np.linalg.norm(self.phys['pore.charge_conservation.rate'])
        assert_allclose(actual=y, desired=1.155203e-13, rtol=1e-5)
        self.phys.add_model(propname='pore.charge_conservation',
                            model=model,
                            phase=self.sw,
                            p_alg=self.p,
                            e_alg=[self.eA, self.eB],
                            assumption='poisson_2d')
        y = np.linalg.norm(self.phys['pore.charge_conservation.rate'])
        assert_allclose(actual=y, desired=2.802853e-7, rtol=1e-5)


if __name__ == '__main__':
    t = MultiphysicsNernstPlanckSolverTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('Running test: '+item)
            t.__getattribute__(item)()

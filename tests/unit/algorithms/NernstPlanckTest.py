import openpnm as op
import scipy as sp
from numpy.testing import assert_allclose


class NernstPlanckTest:

    def setup_class(self):
        sp.random.seed(0)
        self.net = op.network.Cubic(shape=[4, 3, 1], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.area'] = 2e-14
        self.geo['throat.area'] = 1e-14
        self.geo['throat.conduit_lengths.pore1'] = 0.1
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.1

        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.diffusivity.ionX'] = 1e-9
        self.phase['throat.valence.ionX'] = 1

        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.diffusive_conductance.ionX'] = 1e-15
        self.phys['throat.hydraulic_conductance'] = 1e-15
        self.phys['throat.ionic_conductance'] = 1e-15

        self.sf = op.algorithms.StokesFlow(network=self.net, phase=self.phase)
        self.sf.set_value_BC(pores=self.net.pores('back'), values=1)
        self.sf.set_value_BC(pores=self.net.pores('front'), values=0)
        self.sf.run()

        self.phase.update(self.sf.results())

        self.p = op.algorithms.OhmicConduction(network=self.net,
                                               phase=self.phase)
        self.p.settings['conductance'] = 'throat.ionic_conductance'
        self.p.settings['quantity'] = 'pore.potential'
        self.p.set_value_BC(pores=self.net.pores('left'), values=0.05)
        self.p.set_value_BC(pores=self.net.pores('right'), values=0.00)
        self.p.run()

        self.phase.update(self.p.results())

        self.adm = op.algorithms.NernstPlanck(network=self.net,
                                              phase=self.phase, ion='ionX')
        self.adm.set_value_BC(pores=self.net.pores('back'), values=2)
        self.adm.set_value_BC(pores=self.net.pores('front'), values=0)

    def test_powerlaw_NernstPlanck(self):
        mod = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(propname='throat.ad_dif_mig_conductance_powerlaw',
                            model=mod, s_scheme='powerlaw', ion='ionX')
        self.phys.regenerate_models()

        self.adm.setup(conductance='throat.ad_dif_mig_conductance_powerlaw')
        self.adm.run()
        x = [0.,      0.,      0.,
             0.89653, 0.89653, 0.89653,
             1.53924, 1.53924, 1.53924,
             2.,      2.,      2.]
        y = sp.around(self.adm['pore.concentration.ionX'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_upwind_NernstPlanck(self):
        mod = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(propname='throat.ad_dif_mig_conductance_upwind',
                            model=mod, s_scheme='upwind', ion='ionX')
        self.phys.regenerate_models()
        self.adm.setup(conductance='throat.ad_dif_mig_conductance_upwind')
        self.adm.run()
        x = [0.,      0.,      0.,
             0.86486, 0.86486, 0.86486,
             1.51351, 1.51351, 1.51351,
             2.,      2.,      2.]
        y = sp.around(self.adm['pore.concentration.ionX'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_hybrid_NernstPlanck(self):
        mod = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(propname='throat.ad_dif_mig_conductance_hybrid',
                            model=mod, s_scheme='hybrid', ion='ionX')
        self.phys.regenerate_models()

        self.adm.setup(conductance='throat.ad_dif_mig_conductance_hybrid')
        self.adm.run()
        x = [0.,      0.,      0.,
             0.89908, 0.89908, 0.89908,
             1.54128, 1.54128, 1.54128,
             2.,      2.,      2.]
        y = sp.around(self.adm['pore.concentration.ionX'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_exponential_NernstPlanck(self):
        mod = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(propname='throat.ad_dif_mig_conductance_exp',
                            model=mod, s_scheme='exponential', ion='ionX')
        self.phys.regenerate_models()

        self.adm.setup(conductance='throat.ad_dif_mig_conductance_exp')
        self.adm.run()
        x = [0.,      0.,      0.,
             0.89688, 0.89688, 0.89688,
             1.53953, 1.53953, 1.53953,
             2.,      2.,      2.]
        y = sp.around(self.adm['pore.concentration.ionX'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':
    t = NernstPlanckTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t

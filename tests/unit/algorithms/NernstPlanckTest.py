import pytest
import numpy as np
import openpnm as op
from numpy.testing import assert_allclose


class NernstPlanckTest:

    def setup_class(self):
        np.random.seed(0)
        self.net = op.network.Cubic(shape=[4, 3, 1], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['throat.conduit_lengths.pore1'] = 0.1
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.1

        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.diffusivity.ionX'] = 1e-9
        self.phase['throat.valence.ionX'] = 1

        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.diffusive_conductance.ionX'] = 1e-16
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
        self.adm.settings.update({"cache_A": False, "cache_b": False})
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
             1.27816, 1.79057, 2.70356,
             1.59724, 1.93331, 2.46112,
             2.,      2.,      2.]
        y = np.around(self.adm['pore.concentration.ionX'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_upwind_NernstPlanck(self):
        mod = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(propname='throat.ad_dif_mig_conductance_upwind',
                            model=mod, s_scheme='upwind', ion='ionX')
        self.phys.regenerate_models()
        self.adm.setup(conductance='throat.ad_dif_mig_conductance_upwind')
        self.adm.run()
        x = [0.,      0.,      0.,
             1.15437, 1.497,   2.02144,
             1.60103, 1.87748, 2.27264,
             2.,      2.,      2.]
        y = np.around(self.adm['pore.concentration.ionX'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_hybrid_NernstPlanck(self):
        mod = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(propname='throat.ad_dif_mig_conductance_hybrid',
                            model=mod, s_scheme='hybrid', ion='ionX')
        self.phys.regenerate_models()

        self.adm.setup(conductance='throat.ad_dif_mig_conductance_hybrid')
        self.adm.run()
        x = [0.,      0.,      0.,
             1.29501, 1.84357, 2.86142,
             1.58876, 1.93152, 2.47971,
             2.,      2.,      2.]
        y = np.around(self.adm['pore.concentration.ionX'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_exponential_NernstPlanck(self):
        mod = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(propname='throat.ad_dif_mig_conductance_exp',
                            model=mod, s_scheme='exponential', ion='ionX')
        self.phys.regenerate_models()

        self.adm.setup(conductance='throat.ad_dif_mig_conductance_exp')
        self.adm.run()
        x = [0.,      0.,      0.,
             1.2788,  1.79357, 2.71385,
             1.59656, 1.93321, 2.46286,
             2.,      2.,      2.]
        y = np.around(self.adm['pore.concentration.ionX'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_outflow_BC(self):
        mod = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(propname='throat.ad_dif_mig_conductance',
                            model=mod, s_scheme='powerlaw', ion='ionX')
        self.phys.regenerate_models()
        adm = self.adm
        adm.setup(conductance='throat.ad_dif_mig_conductance')
        adm.remove_BC()
        adm.set_value_BC(pores=self.net.pores('back'), values=1)
        adm.set_outflow_BC(pores=self.net.pores('front'))
        adm.set_value_BC(pores=self.net.pores('left'), values=0.1)
        adm.set_value_BC(pores=self.net.pores('right'), values=0.1)
        adm.run()
        x = [0.06676, 0.36184, 0.05913,
             0.1,     0.43529, 0.1,
             0.1,     0.64825, 0.1,
             0.1,     1.,      0.1]
        y = np.around(self.adm['pore.concentration.ionX'], decimals=5)
        assert_allclose(actual=y, desired=x)

    def test_unsupported_scheme_NernstPlanck(self):
        mod = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        with pytest.raises(Exception):
            self.phys.add_model(propname='throat.ad_dif_mig_conductance_exp',
                                model=mod, s_scheme='unsupported_scheme',
                                ion='ionX')

    def test_ad_dif_mig_cond_w_Nt_by_2_dif_cond(self):
        gd = self.phase["throat.diffusive_conductance.ionX"]
        self.phys["throat.Nt_by_2.ionX"] = np.vstack((gd, gd)).T
        mod = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        self.phys.add_model(propname='throat.ad_dif_mig_conductance_Nt_by_2',
                            model=mod, s_scheme='upwind', ion='ionX',
                            throat_diffusive_conductance="throat.Nt_by_2")
        gd_old = self.phys["throat.ad_dif_mig_conductance_upwind"]
        gd_new = self.phys["throat.ad_dif_mig_conductance_Nt_by_2"]
        # New conductance based on (Nt,2) dif_cond must match the old values
        assert_allclose(actual=gd_new, desired=gd_old)

    def test_ad_dif_mig_cond_when_dif_cond_in_wrong_shape(self):
        gd = self.phase["throat.diffusive_conductance.ionX"]
        self.phys["throat.Nt_by_3.ionX"] = np.vstack((gd, gd, gd)).T
        mod = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
        with pytest.raises(Exception):
            self.phys.add_model(
                    propname='throat.ad_dif_mig_conductance_Nt_by_2',
                    model=mod, s_scheme='upwind', ion='ionX',
                    throat_diffusive_conductance="throat.Nt_by_3")

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':
    t = NernstPlanckTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

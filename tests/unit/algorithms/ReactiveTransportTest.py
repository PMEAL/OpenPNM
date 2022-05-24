import pytest
import openpnm as op
import numpy as np
from openpnm.models.physics import source_terms
from numpy.testing import assert_allclose
import copy


class ReactiveTransportTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['throat.diffusive_conductance'] = 1e-15
        self.phase['pore.A'] = -1e-15
        self.phase['pore.k'] = 2
        self.phase.add_model(
            propname='pore.reaction', model=source_terms.standard_kinetics,
            prefactor='pore.A', exponent='pore.k',
            X='pore.concentration', regen_mode='deferred')
        self.alg = op.algorithms.ReactiveTransport(network=self.net, phase=self.phase)

    def test_settings(self):
        temp = copy.deepcopy(self.alg.settings)
        self.alg.settings._update({'conductance': "throat.cond",
                                   'quantity': "pore.test",
                                   'newton_maxiter': 123,
                                   'relaxation_quantity': 3.21})
        assert self.alg.settings["conductance"] == "throat.cond"
        assert self.alg.settings["quantity"] == "pore.test"
        assert self.alg.settings["newton_maxiter"] == 123
        assert self.alg.settings["relaxation_quantity"] == 3.21
        self.alg.settings = temp

    # def test_get_iterative_props(self):
    #     alg = op.algorithms.ReactiveTransport(network=self.net, phase=self.phase)
    #     # When quantity is None
    #     iterative_props = alg._get_iterative_props()
    #     assert len(iterative_props) == 0
    #     # When there's no dependent property
    #     self.alg.settings["quantity"] = "pore.foo"
    #     iterative_props = alg._get_iterative_props()
    #     assert len(iterative_props) == 0
    #     # When there's one dependent property
    #     self.phase.add_model(propname="pore.bar_depends_on_foo",
    #                          model=lambda target, bar="pore.foo": 0.0)
    #     iterative_props = alg._get_iterative_props()
    #     assert len(iterative_props) == 1
    #     assert "pore.bar_depends_on_foo" in iterative_props
    #     # When there are multiple dependent properties (direct and indirect)
    #     self.phase.add_model(propname="pore.baz_depends_on_bar",
    #                          model=lambda target, bar="pore.bar_depends_on_foo": 0.0)
    #     iterative_props = alg._get_iterative_props()
    #     assert len(iterative_props) == 2
    #     assert "pore.baz_depends_on_bar" in iterative_props

    def test_multiple_set_source_with_same_name_should_only_keep_one(self):
        self.alg.settings._update({'conductance': 'throat.diffusive_conductance',
                                   'quantity': 'pore.concentration'})
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        self.alg.run()
        c_mean_desired = 0.717129
        c_mean = self.alg['pore.concentration'].mean()
        assert_allclose(c_mean, c_mean_desired, rtol=1e-5)

    def test_one_value_one_source(self):
        self.alg['pore.bc_rate'] = np.nan
        self.alg['pore.bc_value'] = np.nan
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        self.alg.run()
        c_mean_desired = 0.717129
        c_mean = self.alg['pore.concentration'].mean()
        assert_allclose(c_mean, c_mean_desired, rtol=1e-6)

    def test_source_over_BCs(self):
        self.alg['pore.bc_rate'] = np.nan
        self.alg['pore.bc_value'] = np.nan
        _ = [self.alg.__setitem__(k, False) for k in self.alg.settings.sources]
        self.alg.set_value_BC(pores=self.net.pores('left'), values=1.0)
        self.alg.set_value_BC(pores=self.net.pores('right'), values=0.5)
        with pytest.raises(Exception):
            self.alg.set_source(pores=self.net.pores(["left", "right"]),
                                propname='pore.reaction')

    def test_BCs_over_source(self):
        self.alg['pore.bc_rate'] = np.nan
        self.alg['pore.bc_value'] = np.nan
        _ = [self.alg.__setitem__(k, False) for k in self.alg.settings.sources]
        self.alg.set_source(pores=self.net.pores('left'), propname='pore.reaction')
        with pytest.raises(Exception):
            self.alg.set_value_BC(pores=self.net.pores('left'), values=1.0)

    def test_multiple_source_terms_same_location(self):
        self.alg['pore.bc_rate'] = np.nan
        self.alg['pore.bc_value'] = np.nan
        _ = [self.alg.__setitem__(k, False) for k in self.alg.settings.sources]
        std_kinetics = op.models.physics.generic_source_term.standard_kinetics
        self.phase.add_model(propname='pore.another_reaction', model=std_kinetics,
                             prefactor='pore.A', exponent='pore.k',
                             X='pore.concentration',
                             regen_mode='deferred')
        self.alg.set_source(pores=self.net.pores('left'), propname='pore.reaction')
        self.alg.set_source(pores=self.net.pores('left'),
                            propname='pore.another_reaction')
        self.alg.set_value_BC(pores=self.net.pores('right'), values=1.0)
        self.alg.run()
        cavg = self.alg["pore.concentration"].mean()
        assert_allclose(cavg, 0.666667, rtol=1e-5)

    def test_source_term_is_set_as_iterative_prop(self):
        self.alg['pore.bc_rate'] = np.nan
        self.alg['pore.bc_value'] = np.nan
        _ = [self.alg.__setitem__(k, False) for k in self.alg.settings.sources]
        self.alg.set_source(pores=self.net.pores('left'), propname='pore.reaction')
        assert "pore.reaction" in self.alg._get_iterative_props()

    def test_quantity_relaxation_consistency_w_base_solution(self):
        self.alg['pore.bc_rate'] = np.nan
        self.alg['pore.bc_value'] = np.nan
        _ = [self.alg.__setitem__(k, False) for k in self.alg.settings.sources]
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        self.alg.run()
        c_mean_base = self.alg['pore.concentration'].mean()
        self.alg.settings['relaxation_quantity'] = 0.5
        self.alg.run()
        c_mean_relaxed = self.alg['pore.concentration'].mean()
        assert_allclose(c_mean_base, c_mean_relaxed, rtol=1e-6)

    def test_set_source_with_modes(self):
        self.alg['pore.bc_rate'] = np.nan
        self.alg['pore.bc_value'] = np.nan
        _ = [self.alg.__setitem__(k, False) for k in self.alg.settings.sources]
        self.alg.set_source(pores=self.net.pores('left'),
                            propname='pore.reaction',
                            mode='overwrite')
        assert self.alg['pore.reaction'].sum() == self.net.num_pores('left')
        self.alg.set_source(pores=[0, 1, 2],
                            propname='pore.reaction',
                            mode='overwrite')
        assert self.alg['pore.reaction'].sum() == 3
        self.alg.set_source(pores=[2, 3, 4], propname='pore.reaction', mode='add')
        assert self.alg['pore.reaction'].sum() == 5

    # def test_remove_source(self):
    #     self.alg['pore.bc_rate'] = np.nan
    #     self.alg['pore.bc_value'] = np.nan
    #     _ = [self.alg.__setitem__(k, False) for k in self.alg.settings.sources]
    #     self.alg.set_source(pores=[2, 3, 4],
    #                         propname='pore.reaction',
    #                         mode='overwrite')
    #     assert self.alg['pore.reaction'].sum() == 3
    #     self.alg.remove_source(propname='pore.reaction', pores=[0, 1])
    #     assert self.alg['pore.reaction'].sum() == 3
    #     self.alg.remove_source(propname='pore.reaction', pores=[0, 2])
    #     assert self.alg['pore.reaction'].sum() == 2
    #     self.alg.remove_source(propname='pore.reaction')
    #     assert 'pore.reaction' not in self.alg.keys()

    def test_source_relaxation_consistency_w_base_solution(self):
        self.alg['pore.bc_rate'] = np.nan
        self.alg['pore.bc_value'] = np.nan
        _ = [self.alg.__setitem__(k, False) for k in self.alg.settings.sources]
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        self.alg.run()
        c_mean_base = self.alg['pore.concentration'].mean()
        self.alg.settings['relaxation_source'] = 0.1
        self.alg.run()
        c_mean_relaxed = self.alg['pore.concentration'].mean()
        assert_allclose(c_mean_base, c_mean_relaxed, rtol=1e-6)

    def test_solution_should_diverge_w_large_relaxation(self):
        self.alg['pore.bc_rate'] = np.nan
        self.alg['pore.bc_value'] = np.nan
        _ = [self.alg.__setitem__(k, False) for k in self.alg.settings.sources]
        self.alg.settings._update({'conductance': 'throat.diffusive_conductance',
                                   'quantity': 'pore.concentration',
                                   'newton_maxiter': 50})
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        self.alg.settings._update({'relaxation_factor': 20.0,
                                   'newton_maxiter': 25})
        self.alg.run()
        assert not self.alg.soln.is_converged

    # FIXME: we no longer want to throw exception when maxiter is reached
    def test_check_divergence_if_maxiter_reached(self):
        self.alg['pore.bc_rate'] = np.nan
        self.alg['pore.bc_value'] = np.nan
        _ = [self.alg.__setitem__(k, False) for k in self.alg.settings.sources]
        self.alg.settings._update({'conductance': 'throat.diffusive_conductance',
                                   'quantity': 'pore.concentration',
                                   'newton_maxiter': 2})
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        self.alg.settings._update({'relaxation_quantity': 1.0,
                                   'newton_maxiter': 2})
        with pytest.raises(Exception):
            raise Exception
        self.alg.settings['newton_maxiter'] = 5000

    # def test_variable_conductance(self):
    #     self.alg.reset(bcs=True, source_terms=True)

    #     # Define concentration-dependent diffusivity
    #     def D_var(target, pore_concentration="pore.concentration"):
    #         X = target[pore_concentration]
    #         return 1e-9 * (1 + 5 * (X / 10) ** 2)

    #     # Define a custom diffusive_conductance model dependent on diffusivity
    #     def g_var(target, pore_diffusivity="pore.diffusivity"):
    #         Dt = target["throat.diffusivity"]
    #         return 1e-6 * Dt

    #     self.alg["pore.concentration"] = 0.0
    #     self.phase.add_model(propname="pore.diffusivity", model=D_var)
    #     self.phase.add_model(propname="throat.diffusive_conductance", model=g_var)
    #     self.alg.set_value_BC(pores=self.net.pores("front"), values=10.0)
    #     self.alg.set_value_BC(pores=self.net.pores("back"), values=0.0)
    #     self.alg.run()
    #     shape = op.topotools.get_shape(self.net)
    #     c_avg = self.alg["pore.concentration"].reshape(shape).mean(axis=(0, 2))
    #     desired = [10.0, 8.18175755, 5.42194391, 0.0]
    #     assert_allclose(c_avg, desired, rtol=1e-5)

    # def test_reset(self):
    #     self.alg['pore.bc_rate'] = np.nan
    #     self.alg['pore.bc_value'] = np.nan
    #     _ = [self.alg.__setitem__(k, False) for k in self.alg.settings.sources]
    #     self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
    #     self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
    #     assert 'sources' in self.alg.settings._attrs
    #     self.alg.reset(source_terms=True)
    #     assert not self.alg.settings['sources']

    # def test_ensure_settings_are_valid(self):
    #     alg = op.algorithms.ReactiveTransport(network=self.net,
    #                                           phase=self.phase)
    #     # Commenting this until settings are settled
    #     with pytest.raises(Exception, match=r".*quantity.*"):
    #         alg.run()
    #     alg.settings['quantity'] = 'pore.concentration'
    #     with pytest.raises(Exception, match=r".*conductance.*"):
    #         alg.run()
    #     alg.settings['conductance'] = 'throat.conductance'
    #     with pytest.raises(Exception):
    #         alg.run()

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = ReactiveTransportTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()

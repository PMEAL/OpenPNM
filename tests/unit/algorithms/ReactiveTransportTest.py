import pytest
import openpnm as op
from numpy.testing import assert_allclose


class ReactiveTransportTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.geo = op.geometry.GenericGeometry(
            network=self.net, pores=self.net.Ps, throats=self.net.Ts
        )
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(
            network=self.net, phase=self.phase, geometry=self.geo
        )
        self.phys['throat.diffusive_conductance'] = 1e-15
        self.phys['pore.A'] = -1e-15
        self.phys['pore.k'] = 2
        std_kinetics = op.models.physics.generic_source_term.standard_kinetics
        self.phys.add_model(
            propname='pore.reaction', model=std_kinetics,
            prefactor='pore.A', exponent='pore.k',
            X='pore.concentration', regen_mode='deferred'
        )
        self.alg = op.algorithms.ReactiveTransport(network=self.net,
                                                   phase=self.phase)

    def test_settings(self):
        temp = self.alg.settings.copy()
        self.alg.settings.update({
            'conductance': "throat.cond",
            'quantity': "pore.test",
            'nlin_max_iter': 123,
            'relaxation_source': 1.23,
            'relaxation_quantity': 3.21
        })
        assert self.alg.settings["conductance"] == "throat.cond"
        assert self.alg.settings["quantity"] == "pore.test"
        assert self.alg.settings["nlin_max_iter"] == 123
        assert self.alg.settings["relaxation_source"] == 1.23
        assert self.alg.settings["relaxation_quantity"] == 3.21
        self.alg.settings = temp

    def test_get_iterative_props(self):
        # When quantity is None
        iterative_props = self.alg._get_iterative_props()
        assert len(iterative_props) == 0
        # When there's no dependent property
        self.alg.settings["quantity"] = "pore.foo"
        iterative_props = self.alg._get_iterative_props()
        assert len(iterative_props) == 0
        # When there's one dependent property
        self.phase.add_model(propname="pore.bar_depends_on_foo",
                             model=lambda target, bar="pore.foo": 0.0)
        iterative_props = self.alg._get_iterative_props()
        assert len(iterative_props) == 1
        assert "pore.bar_depends_on_foo" in iterative_props
        # When there are multiple dependent properties (direct and indirect)
        self.phys.add_model(propname="pore.baz_depends_on_bar",
                            model=lambda target, bar="pore.bar_depends_on_foo": 0.0)
        iterative_props = self.alg._get_iterative_props()
        assert len(iterative_props) == 2
        assert "pore.baz_depends_on_bar" in iterative_props

    def test_multiple_set_source_with_same_name_should_only_keep_one(self):
        self.alg.settings.update({'conductance': 'throat.diffusive_conductance',
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
        self.alg.reset(bcs=True, source_terms=True)
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        self.alg.run()
        c_mean_desired = 0.717129
        c_mean = self.alg['pore.concentration'].mean()
        assert_allclose(c_mean, c_mean_desired, rtol=1e-6)

    def test_source_over_BCs(self):
        self.alg.reset(bcs=True, source_terms=True)
        self.alg.set_value_BC(pores=self.net.pores('left'), values=1.0)
        self.alg.set_value_BC(pores=self.net.pores('right'), values=0.5)
        with pytest.raises(Exception):
            self.alg.set_source(pores=self.net.pores(["left", "right"]),
                                propname='pore.reaction')

    def test_BCs_over_source(self):
        self.alg.reset(bcs=True, source_terms=True)
        self.alg.set_source(pores=self.net.pores('left'), propname='pore.reaction')
        with pytest.raises(Exception):
            self.alg.set_value_BC(pores=self.net.pores('left'), values=1.0)

    def test_multiple_source_terms_same_location(self):
        self.alg.reset(bcs=True, source_terms=True)
        std_kinetics = op.models.physics.generic_source_term.standard_kinetics
        self.phys.add_model(propname='pore.another_reaction', model=std_kinetics,
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
        self.alg.reset(bcs=True, source_terms=True)
        self.alg.set_source(pores=self.net.pores('left'), propname='pore.reaction')
        assert "pore.reaction" in self.alg._get_iterative_props()

    def test_quantity_relaxation_consistency_w_base_solution(self):
        self.alg.reset(bcs=True, source_terms=True)
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        self.alg.run()
        c_mean_base = self.alg['pore.concentration'].mean()
        self.alg.settings['relaxation_quantity'] = 0.5
        self.alg.run()
        c_mean_relaxed = self.alg['pore.concentration'].mean()
        assert_allclose(c_mean_base, c_mean_relaxed, rtol=1e-6)

    def test_set_source_with_modes(self):
        self.alg.reset(bcs=True, source_terms=True)
        self.alg.set_source(pores=self.net.pores('left'),
                            propname='pore.reaction',
                            mode='overwrite')
        assert self.alg['pore.reaction'].sum() == self.net.num_pores('left')
        self.alg.set_source(pores=[0, 1, 2],
                            propname='pore.reaction',
                            mode='overwrite')
        assert self.alg['pore.reaction'].sum() == 3
        self.alg.set_source(pores=[2, 3, 4],
                            propname='pore.reaction',
                            mode='merge')
        assert self.alg['pore.reaction'].sum() == 5

    def test_remove_source(self):
        self.alg.reset(bcs=True, source_terms=True)
        self.alg.set_source(pores=[2, 3, 4],
                            propname='pore.reaction',
                            mode='overwrite')
        assert self.alg['pore.reaction'].sum() == 3
        self.alg.remove_source(propname='pore.reaction', pores=[0, 1])
        assert self.alg['pore.reaction'].sum() == 3
        self.alg.remove_source(propname='pore.reaction', pores=[0, 2])
        assert self.alg['pore.reaction'].sum() == 2
        self.alg.remove_source(propname='pore.reaction')
        assert self.alg['pore.reaction'].sum() == 0

    def test_source_relaxation_consistency_w_base_solution(self):
        self.alg.reset(bcs=True, source_terms=True)
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        self.alg.run()
        c_mean_base = self.alg['pore.concentration'].mean()
        self.alg.settings['relaxation_source'] = 0.1
        self.alg.run()
        c_mean_relaxed = self.alg['pore.concentration'].mean()
        assert_allclose(c_mean_base, c_mean_relaxed, rtol=1e-6)

    def test_solution_should_diverge_w_large_relaxation(self):
        self.alg.reset(bcs=True, source_terms=True)
        self.alg.settings.update({'conductance': 'throat.diffusive_conductance',
                                  'quantity': 'pore.concentration',
                                  'nlin_max_iter': 25})
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        self.alg.settings['relaxation_quantity'] = 20.0
        self.alg.settings['relaxation_source'] = 1.0
        with pytest.raises(Exception):
            self.alg.run()
        self.alg.settings['relaxation_quantity'] = 1.0
        self.alg.settings['relaxation_source'] = 20.0
        with pytest.raises(Exception):
            self.alg.run()

    def test_check_divergence_if_maxiter_reached(self):
        self.alg.reset(bcs=True, source_terms=True)
        self.alg.settings.update({'conductance': 'throat.diffusive_conductance',
                                  'quantity': 'pore.concentration',
                                  'nlin_max_iter': 2})
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        self.alg.settings['relaxation_quantity'] = 1.0
        self.alg.settings['relaxation_source'] = 1.0
        with pytest.raises(Exception):
            self.alg.run()
        self.alg.settings['nlin_max_iter'] = 5000

    def test_variable_conductance(self):
        self.alg.reset(bcs=True, source_terms=True)

        # Define concentration-dependent diffusivity
        def variable_diffusivity(target, pore_concentration="pore.concentration"):
            X = target[pore_concentration]
            return 1e-9 * (1 + 5 * (X / 10) ** 2)

        # Define a custom diffusive_conductance model dependent on diffusivity
        def variable_conductance(target, pore_diffusivity="pore.diffusivity"):
            Dt = target["throat.diffusivity"]
            return 1e-6 * Dt

        self.alg["pore.concentration"] = 0.0
        self.phase.add_model(propname="pore.diffusivity",
                             model=variable_diffusivity)
        self.phys.add_model(propname="throat.diffusive_conductance",
                            model=variable_conductance)
        self.alg.set_value_BC(pores=self.net.pores("front"), values=10.0)
        self.alg.set_value_BC(pores=self.net.pores("back"), values=0.0)
        self.alg.run()
        shape = op.topotools.get_shape(self.net)
        c_avg = self.alg["pore.concentration"].reshape(shape).mean(axis=(0, 2))
        desired = [10.0, 8.18175755, 5.42194391, 0.0]
        assert_allclose(c_avg, desired)

    def test_reset(self):
        self.alg.reset(bcs=True, source_terms=True)
        self.alg.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        self.alg.set_value_BC(pores=self.net.pores('top'), values=1.0)
        assert 'sources' in self.alg.settings.keys()
        self.alg.reset(source_terms=True)
        assert not self.alg.settings['sources']

    def test_ensure_settings_are_valid(self):
        alg = op.algorithms.ReactiveTransport(network=self.net,
                                              phase=self.phase)
        with pytest.raises(Exception, match=r".*quantity.*"):
            alg.run()
        alg.settings['quantity'] = 'pore.concentration'
        with pytest.raises(Exception, match=r".*conductance.*"):
            alg.run()
        alg.settings['conductance'] = 'throat.conductance'
        with pytest.raises(Exception):
            alg.run()

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = ReactiveTransportTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: ' + item)
            t.__getattribute__(item)()
    self = t

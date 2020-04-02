import pytest
import openpnm as op
from numpy.testing import assert_allclose


class ReactiveTransportTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[9, 9, 9])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.diffusive_conductance'] = 1e-15
        self.phys['pore.A'] = -1e-15
        self.phys['pore.k'] = 2
        std_kinetics = op.models.physics.generic_source_term.standard_kinetics
        self.phys.add_model(propname='pore.reaction', model=std_kinetics,
                            prefactor='pore.A', exponent='pore.k',
                            X='pore.concentration',
                            regen_mode='deferred')

    def test_set_variable_props(self):
        alg = op.algorithms.ReactiveTransport(network=self.net,
                                              phase=self.phase)
        assert len(alg.settings["variable_props"]) == 0
        alg.set_variable_props(propnames="pore.pressure")
        assert "pore.pressure" in alg.settings["variable_props"]
        # Ensure each prop is only added once
        alg.set_variable_props(propnames="pore.pressure")
        assert len(alg.settings["variable_props"]) == 1
        alg.set_variable_props(propnames=["pore.temperature", "pore.pressure"])
        assert len(alg.settings["variable_props"]) == 2
        assert "pore.pressure" in alg.settings["variable_props"]
        assert "pore.temperature" in alg.settings["variable_props"]

    def test_find_iterative_props(self):
        alg = op.algorithms.ReactiveTransport(network=self.net,
                                              phase=self.phase)
        # When quantity is None
        iterative_props = alg._find_iterative_props()
        assert len(iterative_props) == 0

        # When there's no dependent property
        alg.settings["quantity"] = "pore.foo"
        iterative_props = alg._find_iterative_props()
        assert len(iterative_props) == 0

        # When there's one dependent property
        def mymodel(target, bar="pore.foo"):
            return 0.0
        self.phase.add_model(propname="pore.bar_depends_on_foo", model=mymodel)
        iterative_props = alg._find_iterative_props()
        assert len(iterative_props) == 1
        assert "pore.bar_depends_on_foo" in iterative_props

        # When there are multiple dependent properties (direct and indirect)
        def mymodel2(target, bar="pore.bar_depends_on_foo"):
            return 0.0
        self.phys.add_model(propname="pore.baz_depends_on_bar", model=mymodel2)
        iterative_props = alg._find_iterative_props()
        assert len(iterative_props) == 2
        assert "pore.baz_depends_on_bar" in iterative_props

        # When settings["variable_props"] is not empty
        def mymodel3(target, bar="pore.blah"):
            return 0.0
        self.phys.add_model(propname="pore.lambda_depends_on_blah", model=mymodel3)
        alg.set_variable_props("pore.blah")
        iterative_props = alg._find_iterative_props()
        assert len(iterative_props) == 3
        assert "pore.lambda_depends_on_blah" in iterative_props

    def test_multiple_set_source_with_same_name_should_only_keep_one(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        rt.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        rt.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        rt.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        rt.set_value_BC(pores=self.net.pores('top'), values=1.0)
        rt.run()
        c_mean_desired = 0.648268
        c_mean = rt['pore.concentration'].mean()
        assert_allclose(c_mean, c_mean_desired, rtol=1e-6)

    def test_one_value_one_source(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.setup(solver_tol=1e-10, max_iter=5000,
                 relaxation_source=1.0, relaxation_quantity=1.0)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        rt.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        rt.set_value_BC(pores=self.net.pores('top'), values=1.0)
        rt.run()
        c_mean_desired = 0.648268
        c_mean = rt['pore.concentration'].mean()
        assert_allclose(c_mean, c_mean_desired, rtol=1e-6)

    def test_source_over_BCs(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        rt.set_value_BC(pores=self.net.pores('left'), values=1.0)
        rt.set_value_BC(pores=self.net.pores('right'), values=0.5)
        with pytest.raises(Exception):
            rt.set_source(pores=self.net.pores(["left", "right"]),
                          propname='pore.reaction')

    def test_BCs_over_source(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        rt.set_source(pores=self.net.pores('left'), propname='pore.reaction')
        with pytest.raises(Exception):
            rt.set_value_BC(pores=self.net.pores('left'), values=1.0)

    def test_multiple_source_terms_same_location(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        std_kinetics = op.models.physics.generic_source_term.standard_kinetics
        self.phys.add_model(propname='pore.another_reaction', model=std_kinetics,
                            prefactor='pore.A', exponent='pore.k',
                            X='pore.concentration',
                            regen_mode='deferred')
        rt.set_source(pores=self.net.pores('left'), propname='pore.reaction')
        rt.set_source(pores=self.net.pores('left'), propname='pore.another_reaction')
        rt.set_value_BC(pores=self.net.pores('right'), values=1.0)
        rt.run()
        cavg = rt["pore.concentration"].mean()
        assert_allclose(cavg, 0.61034780)

    def test_source_term_is_set_as_iterative_prop(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        rt.set_source(pores=self.net.pores('left'), propname='pore.reaction')
        assert "pore.reaction" in rt._find_iterative_props()

    def test_quantity_relaxation_consistency_w_base_solution(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.setup(solver_tol=1e-10, max_iter=5000,
                 relaxation_source=1.0, relaxation_quantity=1.0)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        rt.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        rt.set_value_BC(pores=self.net.pores('top'), values=1.0)
        rt.run()
        c_mean_base = rt['pore.concentration'].mean()
        rt.settings['relaxation_quantity'] = 0.5
        rt.run()
        c_mean_relaxed = rt['pore.concentration'].mean()
        assert_allclose(c_mean_base, c_mean_relaxed, rtol=1e-6)

    def test_source_relaxation_consistency_w_base_solution(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.setup(solver_tol=1e-10, max_iter=5000,
                 relaxation_source=1.0, relaxation_quantity=1.0)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        rt.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        rt.set_value_BC(pores=self.net.pores('top'), values=1.0)
        rt.run()
        c_mean_base = rt['pore.concentration'].mean()
        rt.settings['relaxation_source'] = 0.1
        rt.run()
        c_mean_relaxed = rt['pore.concentration'].mean()
        assert_allclose(c_mean_base, c_mean_relaxed, rtol=1e-6)

    def test_solution_should_diverge_w_large_relaxation(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.setup(solver_tol=1e-10, max_iter=100,
                 relaxation_source=1.0, relaxation_quantity=1.0)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        rt.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        rt.set_value_BC(pores=self.net.pores('top'), values=1.0)
        rt.settings['relaxation_quantity'] = 2.0
        rt.settings['relaxation_source'] = 1.0
        with pytest.raises(Exception):
            rt.run()
        rt.settings['relaxation_quantity'] = 1.0
        rt.settings['relaxation_source'] = 2.0
        with pytest.raises(Exception):
            rt.run()

    def test_reset(self):
        rt = op.algorithms.ReactiveTransport(network=self.net,
                                             phase=self.phase)
        rt.setup(solver_tol=1e-10, max_iter=5000,
                 relaxation_source=1.0, relaxation_quantity=1.0)
        rt.settings.update({'conductance': 'throat.diffusive_conductance',
                            'quantity': 'pore.concentration'})
        rt.set_source(pores=self.net.pores('bottom'), propname='pore.reaction')
        rt.set_value_BC(pores=self.net.pores('top'), values=1.0)
        assert 'sources' in rt.settings.keys()
        rt.reset(source_terms=True)
        assert not rt.settings['sources']
        rt.set_variable_props(["pore.blah", "pore.foo"])
        assert len(rt.settings["variable_props"]) == 2
        rt.reset(variable_props=True)
        assert not rt.settings["variable_props"]

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = ReactiveTransportTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

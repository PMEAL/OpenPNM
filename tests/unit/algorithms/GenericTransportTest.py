import pytest
import openpnm as op
import scipy as sp


class GenericTransportTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[9, 9, 9])
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)
        self.phase = op.phases.Air(network=self.net)
        self.phase['pore.mole_fraction'] = 0
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.diffusive_conductance'] = 1.0

    def test_results(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        with pytest.raises(Exception):
            alg.results()

    def test_set_solver(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        # Store old values
        family = alg.settings["solver_family"]
        stype = alg.settings["solver_type"]
        tol = alg.settings["solver_tol"]
        atol = alg.settings["solver_atol"]
        rtol = alg.settings["solver_rtol"]
        maxiter = alg.settings["solver_maxiter"]
        # Set solver settings, but don't provide any arguments
        alg.set_solver()
        # Make sure nothing was changed
        assert alg.settings["solver_family"] == family
        assert alg.settings["solver_type"] == stype
        assert alg.settings["solver_tol"] == tol
        assert alg.settings["solver_atol"] == atol
        assert alg.settings["solver_rtol"] == rtol
        assert alg.settings["solver_maxiter"] == maxiter
        # Set solver settings, this time change everything
        alg.set_solver(solver_family="petsc", solver_type="gmres", maxiter=13,
                       preconditioner="ilu", tol=1e-3, atol=1e-12, rtol=1e-2)
        # Make changes went through
        assert alg.settings["solver_family"] == "petsc"
        assert alg.settings["solver_type"] == "gmres"
        assert alg.settings["solver_tol"] == 1e-3
        assert alg.settings["solver_atol"] == 1e-12
        assert alg.settings["solver_rtol"] == 1e-2
        assert alg.settings["solver_maxiter"] == 13

    def test_remove_boundary_conditions(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.set_value_BC(pores=self.net.pores('top'), values=1)
        alg.set_value_BC(pores=self.net.pores('bottom'), values=0)
        assert sp.sum(sp.isfinite(alg['pore.bc_value'])) > 0
        alg.remove_BC(pores=self.net.pores('top'))
        assert sp.sum(sp.isfinite(alg['pore.bc_value'])) > 0
        alg.remove_BC(pores=self.net.pores('bottom'))
        assert sp.sum(sp.isfinite(alg['pore.bc_value'])) == 0

    def test_generic_transport(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_value_BC(pores=self.net.pores('top'), values=1)
        alg.set_value_BC(pores=self.net.pores('bottom'), values=0)
        alg.run()

    def test_pestc_wrapper(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_value_BC(pores=self.net.pores('top'), values=1)
        alg.set_value_BC(pores=self.net.pores('bottom'), values=0)
        x = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
        # Test different solvers
        solver_types = ['mumps', 'cg', 'gmres', 'bicg']
        # PETSc is not by default installed, so testing should be optional too
        try:
            import petsc4py
            for solver_type in solver_types:
                alg.set_solver(solver_family="petsc", solver_type=solver_type)
                alg.run()
                y = sp.unique(sp.around(alg['pore.mole_fraction'], decimals=3))
                assert sp.all(x == y)
        except ModuleNotFoundError:
            pass

    def test_two_value_conditions(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_value_BC(pores=self.net.pores('top'), values=1)
        alg.set_value_BC(pores=self.net.pores('bottom'), values=0)
        alg.run()
        x = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
        y = sp.unique(sp.around(alg['pore.mole_fraction'], decimals=3))
        assert sp.all(x == y)

    def test_two_value_conditions_cg(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_value_BC(pores=self.net.pores('top'), values=1)
        alg.set_value_BC(pores=self.net.pores('bottom'), values=0)
        alg.settings['solver'] = 'cg'
        alg.run()
        x = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
        y = sp.unique(sp.around(alg['pore.mole_fraction'], decimals=3))
        assert sp.all(x == y)

    def test_one_value_one_rate(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=self.net.pores('bottom'), values=1)
        alg.set_value_BC(pores=self.net.pores('top'), values=0)
        alg.run()
        x = [0., 1., 2., 3., 4., 5., 6., 7., 8.]
        y = sp.unique(sp.around(alg['pore.mole_fraction'], decimals=3))
        assert sp.all(x == y)

    def test_set_iterative_props(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                              phase=self.phase)
        assert len(alg.settings["iterative_props"]) == 0
        alg.set_iterative_props(propnames="pore.pressure")
        assert "pore.pressure" in alg.settings["iterative_props"]
        # Ensure each prop is only added once
        alg.set_iterative_props(propnames="pore.pressure")
        assert len(alg.settings["iterative_props"]) == 1
        alg.set_iterative_props(propnames=["pore.temperature", "pore.pressure"])
        assert len(alg.settings["iterative_props"]) == 2
        assert "pore.pressure" in alg.settings["iterative_props"]
        assert "pore.temperature" in alg.settings["iterative_props"]

    def test_cache_A(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=self.net.pores('bottom'), values=1)
        alg.set_value_BC(pores=self.net.pores('top'), values=0)
        alg.settings["cache_A"] = True
        alg._build_A()
        x = alg._A.mean()
        self.phys["throat.diffusive_conductance"][1] = 50.0
        alg._build_A()
        y = alg._A.mean()
        # When cache_A is True, A is not recomputed, hence x == y
        assert x == y
        alg.settings["cache_A"] = False
        alg._build_A()
        y = alg._A.mean()
        # When cache_A is False, A must be recomputed, hence x!= y
        assert x != y
        # Revert back changes to objects
        self.setup_class()

    def test_rate_single(self):
        alg = op.algorithms.ReactiveTransport(network=self.net,
                                              phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=self.net.pores("left"), values=1.235)
        alg.set_value_BC(pores=self.net.pores("right"), values=0.0)
        alg.run()
        rate = alg.rate(pores=self.net.pores("right"))[0]
        assert sp.isclose(rate, -1.235*self.net.pores("right").size)
        # Net rate must always be zero at steady state conditions
        assert sp.isclose(alg.rate(pores=self.net.Ps), 0.0)

    def test_rate_multiple(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=[0, 1, 2, 3], values=1.235)
        # Note that pore = 0 is assigned two rate values (rate = sum(rates))
        alg.set_rate_BC(pores=[5, 6, 19, 35, 0], values=3.455)
        alg.set_value_BC(pores=[50, 51, 52, 53], values=0.0)
        alg.run()
        rate = alg.rate(pores=[50, 51, 52, 53])[0]
        assert sp.isclose(rate, -(1.235*4 + 3.455*5))   # 4, 5 are number of pores
        # Net rate must always be zero at steady state conditions
        assert sp.isclose(alg.rate(pores=self.net.Ps), 0.0)

    def test_rate_Nt_by_2_conductance(self):
        net = op.network.Cubic(shape=[1, 6, 1])
        geom = op.geometry.StickAndBall(network=net)
        air = op.phases.Air(network=net)
        water = op.phases.Water(network=net)
        m = op.phases.MultiPhase(phases=[air, water], project=net.project)
        m.set_occupancy(phase=air, pores=[0, 1, 2])
        m.set_occupancy(phase=water, pores=[3, 4, 5])
        const = op.models.misc.constant
        K_water_air = 0.5
        m.set_binary_partition_coef(propname="throat.partition_coef",
                                    phases=[water, air],
                                    model=const, value=K_water_air)
        m._set_automatic_throat_occupancy()
        _ = op.physics.Standard(network=net, phase=m, geometry=geom)
        alg = op.algorithms.GenericTransport(network=net, phase=m)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=0, values=1.235)
        alg.set_value_BC(pores=5, values=0.0)
        alg.run()
        rate = alg.rate(pores=5)[0]
        assert sp.isclose(rate, -1.235)
        # Rate at air-water interface throat (#2) must match imposed rate
        rate = alg.rate(throats=2)[0]
        assert sp.isclose(rate, 1.235)
        # Rate at interface pores (#2 @ air-side, #3 @ water-side) must be 0
        rate_air_side = alg.rate(pores=2)[0]
        rate_water_side = alg.rate(pores=3)[0]
        assert sp.isclose(rate_air_side, 0.0)
        assert sp.isclose(rate_water_side, 0.0)
        # Net rate must always be zero at steady state conditions
        assert sp.isclose(alg.rate(pores=net.Ps), 0.0)

    def test_reset(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=self.net.pores('bottom'), values=1)
        alg.set_value_BC(pores=self.net.pores('top'), values=0)
        alg.run()
        assert ~sp.all(sp.isnan(alg['pore.bc_value']))
        assert ~sp.all(sp.isnan(alg['pore.bc_rate']))
        assert 'pore.mole_fraction' in alg.keys()
        alg.reset(bcs=True, results=False)
        assert sp.all(sp.isnan(alg['pore.bc_value']))
        assert sp.all(sp.isnan(alg['pore.bc_rate']))
        assert 'pore.mole_fraction' in alg.keys()
        alg.reset(bcs=True, results=True)
        assert 'pore.mole_fraction' not in alg.keys()
        alg.set_rate_BC(pores=self.net.pores('bottom'), values=1)
        alg.set_value_BC(pores=self.net.pores('top'), values=0)
        alg.run()

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = GenericTransportTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t

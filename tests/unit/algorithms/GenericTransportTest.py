import pytest
import numpy as np
import openpnm as op
import numpy.testing as nt


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
        max_iter = alg.settings["solver_max_iter"]
        # Set solver settings, but don't provide any arguments
        alg.set_solver()
        # Make sure nothing was changed
        assert alg.settings["solver_family"] == family
        assert alg.settings["solver_type"] == stype
        assert alg.settings["solver_tol"] == tol
        assert alg.settings["solver_atol"] == atol
        assert alg.settings["solver_rtol"] == rtol
        assert alg.settings["solver_max_iter"] == max_iter
        # Set solver settings, this time change everything
        alg.set_solver(solver_family="petsc", solver_type="gmres", max_iter=13,
                       preconditioner="ilu", tol=1e-3, atol=1e-12, rtol=1e-2)
        # Make changes went through
        assert alg.settings["solver_family"] == "petsc"
        assert alg.settings["solver_type"] == "gmres"
        assert alg.settings["solver_tol"] == 1e-3
        assert alg.settings["solver_atol"] == 1e-12
        assert alg.settings["solver_rtol"] == 1e-2
        assert alg.settings["solver_max_iter"] == 13

    def test_remove_boundary_conditions(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.set_value_BC(pores=self.net.pores('top'), values=1)
        alg.set_value_BC(pores=self.net.pores('bottom'), values=0)
        assert np.sum(np.isfinite(alg['pore.bc_value'])) > 0
        alg.remove_BC(pores=self.net.pores('top'))
        assert np.sum(np.isfinite(alg['pore.bc_value'])) > 0
        alg.remove_BC(pores=self.net.pores('bottom'))
        assert np.sum(np.isfinite(alg['pore.bc_value'])) == 0

    def test_generic_transport(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_value_BC(pores=self.net.pores('top'), values=1)
        alg.set_value_BC(pores=self.net.pores('bottom'), values=0)
        alg.run()

    def test_two_value_conditions(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_value_BC(pores=self.net.pores('top'), values=1)
        alg.set_value_BC(pores=self.net.pores('bottom'), values=0)
        alg.run()
        x = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
        y = np.unique(np.around(alg['pore.mole_fraction'], decimals=3))
        assert np.all(x == y)

    def test_one_value_one_rate(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=self.net.pores('bottom'), values=1)
        alg.set_value_BC(pores=self.net.pores('top'), values=0)
        alg.run()
        x = [0., 1., 2., 3., 4., 5., 6., 7., 8.]
        y = np.unique(np.around(alg['pore.mole_fraction'], decimals=3))
        assert np.all(x == y)

    def test_cache_A(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=self.net.pores('bottom'), values=1)
        alg.set_value_BC(pores=self.net.pores('top'), values=0)
        alg.settings["cache_A"] = True
        alg.run()
        x_before = alg["pore.mole_fraction"].mean()
        self.phys["throat.diffusive_conductance"][1] = 50.0
        alg.run()
        x_after = alg["pore.mole_fraction"].mean()
        # When cache_A is True, A is not recomputed, hence x == y
        assert x_before == x_after
        alg.settings["cache_A"] = False
        alg.run()
        x_after = alg["pore.mole_fraction"].mean()
        # When cache_A is False, A must be recomputed, hence x!= y
        assert x_before != x_after
        # Revert back changes to objects
        self.setup_class()

    def test_rate_single_pore(self):
        alg = op.algorithms.ReactiveTransport(network=self.net,
                                              phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        pores = self.net.pores("left")
        alg.set_rate_BC(pores=pores, values=1.235*np.ones(pores.size))
        alg.set_value_BC(pores=self.net.pores("right"), values=0.0)
        alg.run()
        rate = alg.rate(pores=self.net.pores("right"))[0]
        assert np.isclose(rate, -1.235*self.net.pores("right").size)
        # Net rate must always be zero at steady state conditions
        assert np.isclose(alg.rate(pores=self.net.Ps), 0.0)

    def test_rate_multiple_pores(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=[0, 1, 2, 3], values=1.235)
        alg.set_rate_BC(pores=[5, 6, 19, 35, 0], values=3.455)
        # Pore 0 is assigned two rate BCs, only the most recent will be kept
        alg.set_value_BC(pores=[50, 51, 52, 53], values=0.0)
        alg.run()
        rate = alg.rate(pores=[50, 51, 52, 53])[0]
        # 3 and 5 are number of pores in each rate BC
        assert np.isclose(rate, -(1.235*3 + 3.455*5))
        # Net rate must always be zero at steady state conditions
        assert np.isclose(alg.rate(pores=self.net.Ps), 0.0)

    def test_rate_multiple_values(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=[0, 1, 2, 3], values=[0, 3.5, 0.4, -12])
        alg.set_value_BC(pores=[50, 51, 52, 53], values=0.0)
        alg.run()
        rate_individual = alg.rate(pores=[0, 1, 2, 3], mode='single')
        rate_net = alg.rate(pores=[0, 1, 2, 3], mode='group')[0]
        nt.assert_allclose(rate_individual, [0, 3.5, 0.4, -12], atol=1e-10)
        nt.assert_allclose(rate_net, sum([0, 3.5, 0.4, -12]))

    def test_rate_Nt_by_2_conductance(self):
        net = op.network.Cubic(shape=[1, 6, 1])
        geom = op.geometry.StickAndBall(network=net,
                                        pores=net.Ps,
                                        throats=net.Ts)
        air = op.phases.Air(network=net)
        water = op.phases.Water(network=net)
        m = op.phases.MultiPhase(phases=[air, water], project=net.project)
        m.set_occupancy(phase=air, pores=[0, 1, 2])
        m.set_occupancy(phase=water, pores=[3, 4, 5])
        const = op.models.misc.constant
        K_water_air = 0.5
        m.set_binary_partition_coef(
            phases=[water, air], model=const, value=K_water_air
        )
        m._set_automatic_throat_occupancy()
        _ = op.physics.Standard(network=net, phase=m, geometry=geom)
        alg = op.algorithms.GenericTransport(network=net, phase=m)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=0, values=1.235)
        alg.set_value_BC(pores=5, values=0.0)
        alg.run()
        rate = alg.rate(pores=5)[0]
        assert np.isclose(rate, -1.235)
        # Rate at air-water interface throat (#2) must match imposed rate
        rate = alg.rate(throats=2)[0]
        assert np.isclose(rate, 1.235)
        # Rate at interface pores (#2 @ air-side, #3 @ water-side) must be 0
        rate_air_side = alg.rate(pores=2)[0]
        rate_water_side = alg.rate(pores=3)[0]
        assert np.isclose(rate_air_side, 0.0)
        assert np.isclose(rate_water_side, 0.0)
        # Net rate must always be zero at steady state conditions
        assert np.isclose(alg.rate(pores=net.Ps), 0.0)

    def test_reset_settings_and_data(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=self.net.pores('bottom'), values=1)
        alg.set_value_BC(pores=self.net.pores('top'), values=0)
        alg.run()
        assert ~np.all(np.isnan(alg['pore.bc_value']))
        assert ~np.all(np.isnan(alg['pore.bc_rate']))
        assert 'pore.mole_fraction' in alg.keys()
        alg.reset(bcs=True, results=False)
        assert np.all(np.isnan(alg['pore.bc_value']))
        assert np.all(np.isnan(alg['pore.bc_rate']))
        assert 'pore.mole_fraction' in alg.keys()
        alg.reset(bcs=True, results=True)
        assert 'pore.mole_fraction' not in alg.keys()
        alg.set_rate_BC(pores=self.net.pores('bottom'), values=1)
        alg.set_value_BC(pores=self.net.pores('top'), values=0)
        alg.run()

    def test_reset_actual_results(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance_temp'
        self.phase['throat.diffusive_conductance_temp'] = 1.0
        alg.settings['quantity'] = 'pore.concentration'
        alg.set_value_BC(pores=self.net.pores('bottom'), values=1)
        alg.set_value_BC(pores=self.net.pores('top'), values=0)
        alg.run()
        m1 = alg.rate(pores=self.net.pores('top'))
        m2 = -alg.rate(pores=self.net.pores('bottom'))
        # This should pass because the alg has only run once
        np.testing.assert_allclose(m1, m2)
        # Now adjust conductance values and re-run
        self.phase['throat.diffusive_conductance_temp'][[0, 1, 2]] *= 0.1
        alg.run()
        m1 = alg.rate(pores=self.net.pores('top'))
        m2 = -alg.rate(pores=self.net.pores('bottom'))
        # The mass won't balance, so the same test will fail
        with pytest.raises(AssertionError):
            np.testing.assert_allclose(m1, m2)
        # Now use reset method
        alg.reset()
        alg.run()
        m1 = alg.rate(pores=self.net.pores('top'))
        m2 = -alg.rate(pores=self.net.pores('bottom'))
        # Now this will pass again
        np.testing.assert_allclose(m1, m2)

    def test_validate_data_health(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.concentration'
        alg.settings['cache_A'] = False
        alg.set_value_BC(pores=self.net.pores('top'), values=1)
        alg.set_value_BC(pores=self.net.pores('bottom'), values=0)
        # Check if the method can catch NaNs in data
        self.phys['throat.diffusive_conductance'][0] = np.nan
        with pytest.raises(Exception):
            alg.run()
        mod = op.models.misc.from_neighbor_pores
        self.phase["pore.seed"] = np.nan
        self.phys.add_model(propname="throat.diffusive_conductance", model=mod,
                            prop="pore.seed", ignore_nans=False)
        with pytest.raises(Exception):
            alg.run()
        self.phase["pore.seed"] = 1
        self.phys.regenerate_models(propnames="throat.diffusive_conductance")
        # Check if the method can catch unhealthy topology
        Ts = self.net.find_neighbor_throats(pores=0)
        op.topotools.trim(self.net, throats=Ts)
        with pytest.raises(Exception):
            alg.run()
        # Reset network back to original
        self.setup_class()

    def test_total_rate(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        h = self.net.check_network_health()
        op.topotools.trim(self.net, pores=h['trim_pores'])
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_rate_BC(pores=[0, 1, 2, 3], total_rate=1)
        alg.set_value_BC(pores=[50, 51, 52, 53], values=0.0)
        alg.run()
        rate_individual = alg.rate(pores=[0, 1, 2, 3], mode='single')
        nt.assert_allclose(rate_individual, [0.25, 0.25, 0.25, 0.25],
                           atol=1e-10)
        # test exceptions that come from adding total_rate feature
        with pytest.raises(Exception):
            alg.set_rate_BC(pores=[0, 1, 2, 3],
                            total_rate=[0.25, 0.25, 0.25, 0.25])
        with pytest.raises(Exception):
            alg.set_rate_BC(pores=[0, 1, 2, 3], rates=1, total_rate=1)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = GenericTransportTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

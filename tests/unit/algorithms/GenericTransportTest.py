import openpnm as op
import scipy as sp
import pytest


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

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
        alg.set_dirichlet_BC(pores=self.net.pores('top'), values=1)
        alg.set_dirichlet_BC(pores=self.net.pores('bottom'), values=0)
        assert sp.sum(alg['pore.dirichlet']) > 0
        alg.remove_BC(pores=self.net.pores('top'))
        assert sp.sum(alg['pore.dirichlet']) > 0
        alg.remove_BC(pores=self.net.pores('bottom'))
        assert sp.sum(alg['pore.dirichlet']) == 0

    def test_generic_transport(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_dirichlet_BC(pores=self.net.pores('top'), values=1)
        alg.set_dirichlet_BC(pores=self.net.pores('bottom'), values=0)
        alg.run()

    def test_two_dirichlet_conditions(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_dirichlet_BC(pores=self.net.pores('top'), values=1)
        alg.set_dirichlet_BC(pores=self.net.pores('bottom'), values=0)
        alg.run()
        x = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
        y = sp.unique(sp.around(alg['pore.mole_fraction'], decimals=3))
        assert sp.all(x == y)

    def test_one_dirichlet_one_neumann(self):
        alg = op.algorithms.GenericTransport(network=self.net,
                                             phase=self.phase)
        alg.settings['conductance'] = 'throat.diffusive_conductance'
        alg.settings['quantity'] = 'pore.mole_fraction'
        alg.set_neumann_BC(pores=self.net.pores('bottom'), values=1)
        alg.set_dirichlet_BC(pores=self.net.pores('top'), values=0)
        alg.run()
        x = [0., 1., 2., 3., 4., 5., 6., 7., 8.]
        y = sp.unique(sp.around(alg['pore.mole_fraction'], decimals=3))
        assert sp.all(x == y)

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = GenericTransportTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t

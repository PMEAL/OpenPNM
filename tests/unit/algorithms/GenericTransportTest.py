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
        self.phys['pore.A'] = 1e6
        self.phys['pore.k'] = 2
        mod = op.models.physics.generic_source_term.standard_kinetics
        self.phys.add_model(propname='pore.reaction',
                            model=mod,
                            prefactor='pore.A',
                            exponent='pore.k',
                            quantity='pore.mole_fraction',
                            regen_mode='normal')

    def test_two_dirichlet_conditions(self):
        alg = op.algorithms.FickianDiffusion(network=self.net,
                                             phase=self.phase)
        alg.set_dirchlet_BC(pores=self.net.pores('top'), values=1)
        alg.set_dirchlet_BC(pores=self.net.pores('bottom'), values=0)
        alg.run()
        x = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
        y = sp.unique(sp.around(alg['pore.mole_fraction'], decimals=3))
        assert sp.all(x == y)

    def test_one_dirichlet_one_neumann(self):
        alg = op.algorithms.FickianDiffusion(network=self.net,
                                             phase=self.phase)
        alg.set_neumann_BC(pores=self.net.pores('bottom'), values=1)
        alg.set_dirchlet_BC(pores=self.net.pores('top'), values=0)
        alg.run()
        x = [0., 1., 2., 3., 4., 5., 6., 7., 8.]
        y = sp.unique(sp.around(alg['pore.mole_fraction'], decimals=3))
        assert sp.all(x == y)

    def test_one_dirichlet_one_source_term(self):
        alg = op.algorithms.FickianDiffusion(network=self.net,
                                             phase=self.phase)
        alg.set_dirchlet_BC(pores=self.net.pores('top'), values=1)
        alg.set_source(pores=self.net.pores('bottom'),
                       propname='pore.reaction')
        alg.run()
        x = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]
        y = sp.unique(sp.around(alg['pore.mole_fraction'], decimals=3))
        assert sp.all(x == y)

    def teardown_class(self):
        mgr = op.Base.Workspace()
        mgr.clear()


if __name__ == '__main__':

    t = GenericTransportTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t

import openpnm as op
import scipy as sp
import importlib
import numpy.testing as nt


class SolversTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.geom = op.geometry.StickAndBall(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geom)
        self.phys['throat.conductance'] = sp.linspace(1, 5, num=self.net.Nt)
        self.alg = op.algorithms.GenericTransport(network=self.net)
        self.alg.settings.update(quantity='pore.x',
                                 conductance='throat.conductance')
        self.alg.setup(phase=self.phase)
        self.alg.set_value_BC(pores=self.net.pores('left'), values=1.0)
        self.alg.set_value_BC(pores=self.net.pores('bottom'), values=0.0)

    def test_scipy_direct(self):
        solvers = ['spsolve']
        self.alg.settings['solver_family'] = 'scipy'
        for solver in solvers:
            self.alg.settings['solver_type'] = solver
            self.alg.run()
            xmean = self.alg['pore.x'].mean()
            nt.assert_allclose(actual=xmean, desired=0.5875950426)

    def test_scipy_iterative(self):
        solvers = ['bicg', 'bicgstab', 'cg', 'cgs', 'qmr', 'gcrotmk',
                   'gmres', 'lgmres']
        self.alg.settings['solver_family'] = 'scipy'
        self.alg.settings['solver_rtol'] = 1e-08
        for solver in solvers:
            self.alg.settings['solver_type'] = solver
            self.alg.run()
            xmean = self.alg['pore.x'].mean()
            nt.assert_allclose(actual=xmean, desired=0.587595, rtol=1e-5)

    def test_scipy_iterative_diverge(self):
        solvers = ['bicg', 'bicgstab', 'cg', 'cgs', 'qmr', 'gcrotmk',
                   'gmres', 'lgmres']
        self.alg.settings.update(solver_family='scipy',
                                 solver_maxiter=1)
        for solver in solvers:
            self.alg.settings['solver_type'] = solver
            with nt.assert_raises(Exception):
                self.alg.run()
        self.alg.settings.update(solver_maxiter=100)

    def test_pyamg(self):
        self.alg.settings['solver_family'] = 'pyamg'
        if importlib.util.find_spec('pyamg') is None:
            with nt.assert_raises(Exception):
                self.alg.run()
        else:
            self.alg.run()
            xmean = self.alg['pore.x'].mean()
            nt.assert_allclose(actual=xmean, desired=0.587595, rtol=1e-5)

    def test_petsc(self):
        self.alg.settings['solver_family'] = 'petsc'
        if importlib.util.find_spec('petsc4py') is None:
            with nt.assert_raises(Exception):
                self.alg.run()
        else:
            self.alg.run()
            xmean = self.alg['pore.x'].mean()
            nt.assert_allclose(actual=xmean, desired=0.587595, rtol=1e-5)

    def test_nonsymmetric_algorithms_w_cg_solver_should_throw_error(self):
        air = op.phases.Air(network=self.net)
        phys = op.physics.Standard(network=self.net, phase=air, geometry=self.geom)
        ad = op.algorithms.AdvectionDiffusion(network=self.net, phase=air)
        ad.set_value_BC(pores=self.net.pores("left"), values=1.0)
        ad.set_value_BC(pores=self.net.pores("right"), values=0.1)
        sf = op.algorithms.StokesFlow(network=self.net, phase=air)
        sf.set_value_BC(pores=self.net.pores("left"), values=1.0)
        sf.set_value_BC(pores=self.net.pores("right"), values=0.0)
        sf.run()
        air.update(sf.results())
        phys.regenerate_models()
        ad.settings.update({"cache_A": False, "cache_b": False, "solver_type": "cg"})
        with nt.assert_raises(Exception):
            ad.run()


if __name__ == '__main__':
    t = SolversTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t

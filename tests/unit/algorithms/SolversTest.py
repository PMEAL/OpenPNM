import pytest
import importlib
import numpy as np
import numpy.testing as nt
import openpnm as op


class SolversTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10])
        self.net.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders)
        self.net.regenerate_models()
        self.phase = op.phase.Phase(network=self.net)
        self.phase['throat.conductance'] = np.linspace(1, 5, num=self.net.Nt)
        self.alg = op.algorithms.Transport(network=self.net, phase=self.phase)
        self.alg.settings._update({'quantity': 'pore.x',
                                   'conductance': 'throat.conductance'})
        self.alg.set_value_BC(pores=self.net.pores('front'), values=1.0)
        self.alg.set_value_BC(pores=self.net.pores('bottom'), values=0.0)

    # def test_solver_not_available(self):
    #     self.alg.settings['solver_family'] = 'not_supported_solver'
    #     with pytest.raises(Exception):
    #         self.alg.run()

    # def test_scipy_direct(self):
    #     solvers = ['spsolve']
    #     self.alg.settings['solver_family'] = 'scipy'
    #     for solver in solvers:
    #         self.alg.settings['solver_type'] = solver
    #         self.alg.run()
    #         xmean = self.alg['pore.x'].mean()
    #         nt.assert_allclose(actual=xmean, desired=0.5875950426)

    # def test_scipy_iterative(self):
    #     solvers = ['bicg', 'bicgstab', 'cg', 'cgs', 'qmr', 'gcrotmk',
    #                 'gmres', 'lgmres']
    #     self.alg.settings['solver_family'] = 'scipy'
    #     self.alg.settings['solver_rtol'] = 1e-08
    #     for solver in solvers:
    #         self.alg.settings['solver_type'] = solver
    #         self.alg.run()
    #         xmean = self.alg['pore.x'].mean()
    #         nt.assert_allclose(actual=xmean, desired=0.587595, rtol=1e-5)

    # def test_scipy_iterative_diverge(self):
    #     solvers = ['bicg', 'bicgstab', 'cg', 'cgs', 'qmr', 'gcrotmk',
    #                 'gmres', 'lgmres']
    #     self.alg.settings._update(solver_family='scipy',
    #                               solver_max_iter=1)
    #     for solver in solvers:
    #         self.alg.settings['solver_type'] = solver
    #         with nt.assert_raises(Exception):
    #             self.alg.run()
    #     self.alg.settings._update(solver_max_iter=5000)

    # def test_pyamg_exception_if_not_found(self):
    #     self.alg.settings['solver_family'] = 'pyamg'
    #     if not importlib.util.find_spec("pyamg"):
    #         with pytest.raises(Exception):
    #             self.alg.run()

    # TODO: Uncomment this test once pyamg/pyamg#273 is fixed
    # @catch_module_not_found
    # def test_pyamg(self):
    #     self.alg.settings['solver_family'] = 'pyamg'
    #     self.alg.run()
    #     xmean = self.alg['pore.x'].mean()
    #     nt.assert_allclose(actual=xmean, desired=0.587595, rtol=1e-5)

    # def test_pypardiso_exception_if_not_found(self):
    #     self.alg.settings['solver_family'] = 'pypardiso'
    #     if not importlib.util.find_spec("pypardiso"):
    #         with pytest.raises(Exception):
    #             self.alg.run()

    # @catch_module_not_found
    # def test_pypardiso(self):
    #     self.alg.settings['solver_family'] = 'pypardiso'
    #     self.alg.run()
    #     xmean = self.alg['pore.x'].mean()
    #     nt.assert_allclose(actual=xmean, desired=0.587595, rtol=1e-5)

    # def test_petsc_exception_if_not_found(self):
    #     self.alg.settings['solver_family'] = 'petsc'
    #     self.alg.settings['solver_type'] = 'cg'
    #     if not importlib.util.find_spec("petsc4py"):
    #         with pytest.raises(Exception):
    #             self.alg.run()

    # @catch_module_not_found
    # def test_petsc_mumps(self):
    #     self.alg.settings['solver_family'] = 'petsc'
    #     self.alg.settings['solver_type'] = 'mumps'
    #     self.alg.run()
    #     xmean = self.alg['pore.x'].mean()
    #     nt.assert_allclose(actual=xmean, desired=0.587595, rtol=1e-5)

    # @catch_module_not_found
    # def test_petsc_iterative(self):
    #     self.alg.settings['solver_family'] = 'petsc'
    #     iterative_solvers = [
    #         'cg', 'groppcg', 'pipecg', 'pipecgrr',
    #         'nash', 'stcg', 'gltr', 'fcg', 'pipefcg', 'gmres', 'pipefgmres', 'fgmres',
    #         'lgmres', 'dgmres', 'pgmres', 'tcqmr', 'bcgs', 'ibcgs', 'fbcgs', 'fbcgsr',
    #         'bcgsl', 'pipebcgs', 'cgs', 'tfqmr', 'cr', 'pipecr',
    #         'bicg', 'minres', 'symmlq', 'lcd', 'gcr', 'pipegcr',
    #     ]
    #     for solver in iterative_solvers:
    #         self.alg.settings['solver_type'] = solver
    #         self.alg.run()
    #         xmean = self.alg['pore.x'].mean()
    #         nt.assert_allclose(actual=xmean, desired=0.587595, rtol=1e-5)


if __name__ == '__main__':
    t = SolversTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()

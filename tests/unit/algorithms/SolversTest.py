import numpy as np
import numpy.testing as nt
import openpnm as op


class SolversTest:

    def setup_class(self):
        np.random.seed(0)
        self.net = op.network.Cubic(shape=[3, 4, 5])
        mods = op.models.collections.geometry.spheres_and_cylinders
        self.net.add_model_collection(mods)
        self.net.regenerate_models()
        self.phase = op.phase.Phase(network=self.net)
        self.phase['throat.conductance'] = np.linspace(1, 5, num=self.net.Nt)
        self.alg = op.algorithms.Transport(network=self.net, phase=self.phase)
        self.alg.settings._update({'quantity': 'pore.x',
                                   'conductance': 'throat.conductance'})
        self.alg.set_value_BC(pores=self.net.pores('front'), values=1)
        self.alg.set_value_BC(pores=self.net.pores('bottom'), values=0, mode='overwrite')

    def test_scipy_spsolve(self):
        solver = op.solvers.ScipySpsolve()
        self.alg.run(solver=solver)
        x = self.alg['pore.x']
        nt.assert_allclose(x.mean(), 0.624134, rtol=1e-5)

    def test_pardiso_spsolve(self):
        solver = op.solvers.PardisoSpsolve()
        self.alg.run(solver=solver)
        x = self.alg['pore.x']
        nt.assert_allclose(x.mean(), 0.624134, rtol=1e-5)

    def test_scipy_cg(self):
        solver = op.solvers.ScipyCG()
        self.alg.run(solver=solver)
        x = self.alg['pore.x']
        nt.assert_allclose(x.mean(), 0.624134, rtol=1e-5)


if __name__ == '__main__':
    t = SolversTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()

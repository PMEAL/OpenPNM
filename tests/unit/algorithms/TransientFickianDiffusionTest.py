import numpy as np
from numpy.testing import assert_allclose
import openpnm as op


class TransientFickianDiffusionTest:

    def setup_class(self):
        np.random.seed(0)
        self.net = op.network.Cubic(shape=[4, 3, 1], spacing=1.0)
        self.net.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders()
        )
        self.net.regenerate_models()
        self.net['pore.volume'] = 1e-14
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['pore.molar_density'] = 55500
        self.phase['throat.diffusive_conductance'] = 1e-15
        self.alg = op.algorithms.TransientFickianDiffusion(network=self.net,
                                                           phase=self.phase)
        self.alg.settings._update({'quantity': 'pore.concentration',
                                   'conductance': 'throat.diffusive_conductance'})
        self.alg.set_value_BC(pores=self.net.pores('right'), values=1)
        self.alg.set_value_BC(pores=self.net.pores('left'), values=0)

    def test_transient_fickian_diffusion_intermediate_time(self):
        self.alg.run(x0=0, tspan=(0, 10))
        desired = 0.40803
        actual = self.alg.x.mean()
        assert_allclose(actual, desired, rtol=1e-5)

    def test_transient_fickian_diffusion_steady_state(self):
        self.alg.run(x0=0, tspan=(0, 200))  # pick a relatively large tend
        desired = (1 + 0) / 2               # steady value = avearge of BCs
        actual = self.alg.x.mean()
        assert_allclose(actual, desired, rtol=1e-5)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = TransientFickianDiffusionTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
    self = t

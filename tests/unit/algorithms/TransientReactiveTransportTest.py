import numpy as np
import numpy.testing as nt
import openpnm as op
from scipy.interpolate import interp1d


class TransientReactiveTransportTest:

    def setup_class(self):
        np.random.seed(0)
        self.net = op.network.Cubic(shape=[3, 3, 1], spacing=1e-6)
        self.net.add_model_collection(op.models.collections.geometry.spheres_and_cylinders)
        self.net.regenerate_models()
        self.net['pore.volume'] = 1e-12
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['pore.A'] = -1e-13
        self.phase['pore.k'] = 2
        self.phase['throat.diffusive_conductance'] = 1e-12
        mod = op.models.physics.generic_source_term.standard_kinetics
        self.phase.add_model(propname='pore.reaction',
                             model=mod,
                             prefactor='pore.A',
                             exponent='pore.k',
                             X='pore.concentration',
                             regen_mode='deferred')
        self.settings = {'conductance': 'throat.diffusive_conductance',
                         'quantity': 'pore.concentration'}
        self.alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                            phase=self.phase,
                                                            settings=self.settings)
        self.alg.settings._update({'quantity': 'pore.concentration',
                                   'conductance': 'throat.diffusive_conductance'})
        self.alg.set_value_BC(pores=self.net.pores('front'), values=2)
        self.alg.set_source(propname='pore.reaction', pores=self.net.pores('back'))

    def test_transient_reactive_transport(self):
        self.alg.run(x0=0, tspan=(0, 1))
        desired = 1.13133
        actual = self.alg.x.mean()
        nt.assert_allclose(actual, desired, rtol=1e-5)

    def test_transient_solution(self):
        out = self.alg.run(x0=0, tspan=(0, 1), saveat=0.1)
        # Test datatype
        from openpnm.algorithms._solution import TransientSolution
        quantity = self.alg.settings['quantity']
        assert isinstance(out[quantity], TransientSolution)
        # Ensure solution object is attached to the algorithm
        assert isinstance(self.alg.soln[quantity], TransientSolution)
        # Test shape
        nt.assert_array_equal(self.alg.soln[quantity].shape, (self.alg.Np, 11))
        # Test stored time points
        nt.assert_array_equal(self.alg.soln[quantity].t, np.arange(0, 1.1, 0.1))
        # Ensure solution is callable (i.e., interpolates intermediate times)
        assert hasattr(out[quantity], "__call__")
        # Test solution interpolation
        f = interp1d(out[quantity].t, out[quantity])
        nt.assert_allclose(f(0.05), out[quantity](0.05))
        # Ensure no extrapolation
        with nt.assert_raises(Exception):
            out(1.01)

    def test_consecutive_runs_preserves_solution(self):
        # Integrate from 0 to 0.3 in two steps
        self.alg.run(x0=0, tspan=(0, 0.1))
        out1 = self.alg.run(x0=self.alg.x, tspan=(0.1, 0.3))
        # Integrate from 0 to 0.3 in one go
        out2 = self.alg.run(x0=0, tspan=(0, 0.3))
        # Ensure the results match
        quantity = self.alg.settings['quantity']
        nt.assert_allclose(out1[quantity](0.3), out2[quantity](0.3), rtol=1e-5)


    def test_adding_bc_over_sources(self):
        with nt.assert_raises(Exception):
            self.alg.set_value_BC(pores=self.net.pores("right"), values=0.3)

    def test_adding_sources_over_bc(self):
        with nt.assert_raises(Exception):
            self.alg.set_source(propname='pore.reaction',
                                pores=self.net.pores('left'))

    def test_ensure_settings_are_valid(self):
        alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                       phase=self.phase)
        # Commenting the next bits until settings are settled
        # with nt.assert_raises_regex(Exception, r".*quantity.*"):
        #     alg.run(x0=0, tspan=(0, 1))
        # alg.settings['quantity'] = 'pore.concentration'
        # with nt.assert_raises_regex(Exception, r".*conductance.*"):
        #     alg.run(x0=0, tspan=(0, 1))
        # alg.settings['conductance'] = 'throat.diffusive_conductance'
        # alg.run(x0=0, tspan=(0, 1))

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = TransientReactiveTransportTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()

import openpnm as op
import numpy as np
import numpy.testing as nt


class IntegratorsTest:

    def setup_class(self):
        np.random.seed(10)
        self.net = op.network.Cubic([4, 4, 1], spacing=1e-4)
        geo_mods = op.models.collections.geometry.spheres_and_cylinders.copy()
        self.net.add_model_collection(geo_mods)
        self.net.regenerate_models()
        self.phase = op.phase.Water(network=self.net)
        self.phase['pore.diffusivity'] = 1e-10
        self.phase['pore.concentration'] = 100
        phys_mods = op.models.collections.physics.basic.copy()
        self.phase.add_model_collection(phys_mods)
        self.phase.regenerate_models()
        rxn = op.models.physics.source_terms.standard_kinetics
        self.phase['pore.reaction_constant'] = 1e-14
        self.phase['pore.reaction_order'] = 1
        self.phase.add_model(propname='pore.source',
                             model=rxn,
                             X='pore.concentration',
                             prefactor='pore.reaction_constant',
                             exponent='pore.reaction_order')
        self.alg = op.algorithms.TransientReactiveTransport(network=self.net,
                                                            phase=self.phase)
        self.alg.settings['quantity'] = 'pore.concentration'
        self.alg.settings['conductance'] = 'throat.diffusive_conductance'
        self.alg.set_source(pores=[0], propname='pore.source', mode='add')

    def test_scipy_RK45(self):
        x0 = np.ones(self.net.Np)*100
        tspan = (0, 10)
        saveat = 1.0
        integrator = op.integrators.ScipyRK45()
        self.alg.run(x0, tspan, saveat, integrator)
        x = self.alg.x
        nt.assert_allclose(x.mean(), 110.97095388, rtol=1e-5)

    def test_scipy_BDF(self):
        x0 = np.ones(self.net.Np)*100
        tspan = (0, 10)
        saveat = 1.0
        integrator = op.integrators.ScipyBDF()
        self.alg.run(x0, tspan, saveat, integrator)
        x = self.alg.x
        nt.assert_allclose(x.mean(), 110.97112100, rtol=1e-5)

    def test_scipy_LSODA(self):
        x0 = np.ones(self.net.Np)*100
        tspan = (0, 10)
        saveat = 1.0
        integrator = op.integrators.ScipyLSODA()
        self.alg.run(x0, tspan, saveat, integrator)
        x = self.alg.x
        nt.assert_allclose(x.mean(), 110.97097456, rtol=1e-5)


if __name__ == '__main__':
    t = IntegratorsTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()

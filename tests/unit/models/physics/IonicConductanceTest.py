import numpy as np
import openpnm as op
from numpy.testing import assert_allclose


class IonicConductanceTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.net['pore.diameter'] = 1.12
        self.net['throat.diameter'] = 0.56
        self.net['pore.area'] = 1.0
        self.net['throat.cross_sectional_area'] = 1.0
        self.net['pore.volume'] = 1.0
        self.net['throat.volume'] = 1.0
        self.net['throat.diffusive_size_factors'] = \
            np.ones([self.net.Nt, 3])*(0.123, 0.981, 0.551)
        self.phase = op.phase.Phase(network=self.net)
        self.phase['pore.permittivity'] = 78.0

    def test_generic_ionic_poisson(self):
        mod = op.models.physics.ionic_conductance.poisson
        self.phase.add_model(propname='throat.ionic_conductance', model=mod)
        self.phase.regenerate_models()
        actual = np.mean(self.phase['throat.ionic_conductance'])
        assert_allclose(actual, desired=6.298849e-11, rtol=1e-5)

    # def test_ionic_conductance_electroneutrality(self):
    #     sw = mixtures.SalineWater(network=self.net)
    #     Na = sw.components[f'Na_{sw.name}']
    #     Cl = sw.components[f'Cl_{sw.name}']
    #     sw[f'pore.concentration.{Cl.name}'] = np.linspace(10, 20, 27) * 0.5
    #     sw[f'pore.concentration.{Na.name}'] = np.linspace(10, 20, 27)
    #     electroneutrality = op.models.physics.ionic_conductance.electroneutrality
    #     self.phase.add_model(propname='throat.ionic_conductance',
    #                          model=electroneutrality, ions=[Na.name, Cl.name])
    #     self.phase.regenerate_models()
    #     actual = np.mean(self.phase['throat.ionic_conductance'])
    #     assert_allclose(actual, desired=0.0116345, rtol=1e-5)


if __name__ == '__main__':

    t = IonicConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('Running test: '+item)
            t.__getattribute__(item)()

import numpy as np
import openpnm as op
from numpy.testing import assert_allclose


class ElectricalConductanceTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.net['pore.diameter'] = 1.0
        self.net['throat.diameter'] = 0.5
        self.net['throat.diffusive_size_factors'] = \
            np.ones([self.net.Nt, 3])*(0.123, 0.981, 0.551)
        self.phase = op.phase.Phase(network=self.net)
        self.phase['pore.electrical_conductivity'] = 1.0

    def test_generic_electrical(self):
        mod = op.models.physics.electrical_conductance.generic_electrical
        self.phase.add_model(propname='throat.electrical_conductance', model=mod)
        self.phase.regenerate_models()
        actual = np.mean(self.phase['throat.electrical_conductance'])
        assert_allclose(actual, desired=0.091205, rtol=1e-5)

    def test_series_resistors(self):
        mod = op.models.physics.electrical_conductance.series_resistors
        self.phase.add_model(propname='throat.electrical_conductance', model=mod)
        self.phase.regenerate_models()
        actual = np.mean(self.phase['throat.electrical_conductance'])
        assert_allclose(actual, desired=0.091205, rtol=1e-5)


if __name__ == '__main__':

    t = ElectricalConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

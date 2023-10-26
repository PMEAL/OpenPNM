import numpy as np
from numpy.testing import assert_allclose
import openpnm.models.physics.electrical_conductance as electrical_conductance

import openpnm as op


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

    def test_conductance_shape(self):
        available_models = [
            getattr(electrical_conductance, model_name)
            for model_name in dir(electrical_conductance)
            if callable(getattr(electrical_conductance, model_name))
            ]
        for model in available_models:
            G = model(phase=self.phase)
            assert_allclose(G.shape, (self.net.Nt, 2), rtol=0)

if __name__ == '__main__':

    t = ElectricalConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()

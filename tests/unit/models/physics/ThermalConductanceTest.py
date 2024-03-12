import numpy as np
from numpy.testing import assert_allclose
import openpnm.models.physics.thermal_conductance as thermal_conductance

import openpnm as op


class ThermalConductanceTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.net['throat.diffusive_size_factors'] = \
            np.ones([self.net.Nt, 3])*(0.4, 0.2, 0.3)
        self.phase = op.phase.Phase(network=self.net)
        self.phase['pore.thermal_conductivity'] = 0.5

    def test_generic_thermal(self):
        mod = op.models.physics.thermal_conductance.generic_thermal
        self.phase.add_model(propname='throat.thermal_conductance', model=mod)
        self.phase.regenerate_models()
        actual = self.phase['throat.thermal_conductance'].mean()
        desired = 1 / (1/(0.4*0.5) + 1/(0.2*0.5) + 1/(0.3*0.5))
        assert_allclose(actual, desired)

    def test_series_resistors(self):
        mod = op.models.physics.thermal_conductance.series_resistors
        self.phase.add_model(propname='throat.thermal_conductance', model=mod)
        self.phase.regenerate_models()
        actual = self.phase['throat.thermal_conductance'].mean()
        desired = 1 / (1/(0.4*0.5) + 1/(0.2*0.5) + 1/(0.3*0.5))
        assert_allclose(actual, desired)
    
    def test_conductance_shape(self):
        available_models = [
            getattr(thermal_conductance, model_name)
            for model_name in dir(thermal_conductance)
            if callable(getattr(thermal_conductance, model_name))
            ]
        for model in available_models:
            G = model(phase=self.phase)
            assert_allclose(G.shape, (self.net.Nt, 2), rtol=0)


if __name__ == '__main__':

    t = ThermalConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()

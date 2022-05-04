import openpnm as op
from numpy.testing import assert_allclose


class ThermalConductanceTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.net['throat.diffusive_size_factors'] = {
            "pore1": 0.4, "throat": 0.2, "pore2": 0.3
        }
        self.phase = op.phase.GenericPhase(network=self.net)
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


if __name__ == '__main__':

    t = ThermalConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

import openpnm as op
from numpy.testing import assert_allclose


class ThermalConductanceTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['throat.diffusive_size_factors'] = {
            "pore1": 0.4, "throat": 0.2, "pore2": 0.3
        }
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['pore.thermal_conductivity'] = 0.5
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_generic_thermal(self):
        mod = op.models.physics.thermal_conductance.generic_thermal
        self.phys.add_model(propname='throat.thermal_conductance', model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.thermal_conductance'].mean()
        desired = 1 / (1/(0.4*0.5) + 1/(0.2*0.5) + 1/(0.3*0.5))
        assert_allclose(actual, desired)


    def test_series_resistors(self):
        mod = op.models.physics.thermal_conductance.series_resistors
        self.phys.add_model(propname='throat.thermal_conductance', model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.thermal_conductance'].mean()
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

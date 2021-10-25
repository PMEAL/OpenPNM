import openpnm as op
from numpy.testing import assert_allclose


class ElectricalConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['throat.diffusive_size_factors'] = {
            "pore1": 0.4, "throat": 0.2, "pore2": 0.3
        }
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.electrical_conductivity'] = 0.5
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_electrical_conductance(self):
        mod = op.models.physics.electrical_conductance.series_resistors
        self.phys.add_model(propname='throat.electrical_conductance', model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.electrical_conductance'].mean()
        desired = 1 / (1/(0.4*0.5) + 1/(0.2*0.5) + 1/(0.3*0.5))
        assert_allclose(actual, desired)

    def test_generic_electrical(self):
        mod = op.models.physics.electrical_conductance.generic_electrical
        self.phys.add_model(propname='throat.electrical_conductance', model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.electrical_conductance']
        desired = 1 / (1/(0.4*0.5) + 1/(0.2*0.5) + 1/(0.3*0.5))
        assert_allclose(actual, desired)


if __name__ == '__main__':

    t = ElectricalConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('Running test: '+item)
            t.__getattribute__(item)()

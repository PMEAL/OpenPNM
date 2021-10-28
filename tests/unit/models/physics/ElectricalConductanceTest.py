import numpy as np
import openpnm as op
from numpy.testing import assert_allclose


class ElectricalConductanceTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.diameter'] = 1.0
        self.geo['throat.diameter'] = 0.5
        self.geo['throat.diffusive_size_factors'] = {
            "pore1": 0.123, "throat": 0.981, "pore2": 0.551
        }
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.electrical_conductivity'] = 1.0
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_generic_electrical(self):
        mod = op.models.physics.electrical_conductance.generic_electrical
        self.phys.add_model(propname='throat.electrical_conductance', model=mod)
        self.phys.regenerate_models()
        actual = np.mean(self.phys['throat.electrical_conductance'])
        assert_allclose(actual, desired=0.091205, rtol=1e-5)

    def test_series_resistors(self):
        mod = op.models.physics.electrical_conductance.series_resistors
        self.phys.add_model(propname='throat.electrical_conductance', model=mod)
        self.phys.regenerate_models()
        actual = np.mean(self.phys['throat.electrical_conductance'])
        assert_allclose(actual, desired=0.091205, rtol=1e-5)


if __name__ == '__main__':

    t = ElectricalConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

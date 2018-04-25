import openpnm as op
from numpy.testing import assert_approx_equal


class ElectricalConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_electrical_conductance(self):
        self.phase['pore.conductivity'] = 1
        self.geo['pore.area'] = 1
        self.geo['pore.diameter'] = 1
        self.geo['throat.area'] = 1
        self.geo['throat.length'] = 0.0001
        f = op.models.physics.electrical_conductance.series_resistors
        self.phys.add_model(propname='throat.electrical_conductance',
                            conductivity='pore.conductivity',
                            model=f)
        self.phys.regenerate_models()
        a = 0.99990001
        assert_approx_equal(self.phys['throat.electrical_conductance'].mean(), a)


if __name__ == '__main__':

    t = ElectricalConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

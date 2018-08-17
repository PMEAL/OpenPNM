import openpnm as op
from numpy.testing import assert_approx_equal


class ElectricalConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[4, 4, 4])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.area'] = 1
        self.geo['throat.area'] = 1
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.electrical_conductivity'] = 1
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def test_electrical_conductance(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.15
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.25
        mod = op.models.physics.electrical_conductance.series_resistors
        self.phys.add_model(propname='throat.electrical_conductance',
                            model=mod)
        self.phys.regenerate_models()
        actual = self.phys['throat.electrical_conductance'].mean()
        assert_approx_equal(actual, desired=1.0)


if __name__ == '__main__':

    t = ElectricalConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

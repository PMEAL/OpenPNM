import openpnm as op
from numpy.testing import assert_approx_equal


class HydraulicConductanceTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['pore.area'] = 1.0
        self.geo['throat.area'] = 0.5
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.viscosity'] = 1e-5
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_hagen_poiseuille(self):
        self.geo['throat.conduit_lengths.pore1'] = 0.25
        self.geo['throat.conduit_lengths.throat'] = 0.6
        self.geo['throat.conduit_lengths.pore2'] = 0.15
        mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
        self.phys.add_model(propname='throat.hydraulic_conductance',
                            model=mod)
        actual = self.phys['throat.hydraulic_conductance'].mean()
        assert_approx_equal(actual, desired=1421.0262776)


if __name__ == '__main__':

    t = HydraulicConductanceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

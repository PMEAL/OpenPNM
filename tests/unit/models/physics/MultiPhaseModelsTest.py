import openpnm as op
import scipy as sp
from numpy.testing import assert_approx_equal
import openpnm.models.physics as pm


class MultiPhaseModelsTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        Ps = self.net.Ps
        Ts = self.net.Ts
        self.geo = op.geometry.GenericGeometry(network=self.net, pores=Ps,
                                               throats=Ts)
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.occupancy'] = sp.ones(self.net.Np)
        self.phase['throat.occupancy'] = sp.ones(self.net.Nt)
        self.phase['pore.occupancy'][[6, 7, 19, 25]] = 0
        self.phase['throat.occupancy'][[1, 2, 3]] = 0
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        sp.random.seed(0)
        self.phys['throat.capillary_pressure'] = 7000*sp.rand(self.phys.Nt)
        self.phys['throat.diffusive_conductance'] = 2

    def test_conduit_conductance_strict(self):
        self.phys.add_model(propname='throat.conduit_conductance',
                            throat_conductance='throat.diffusive_conductance',
                            model=pm.multiphase.conduit_conductance,
                            mode='strict', factor=0)
        Tinv = [1, 2, 3, 4, 5, 12, 13, 16, 17, 21, 22, 31, 34, 42, 43, 46, 52]
        a = sp.where(self.phase['throat.conduit_conductance'] == 0)[0]
        assert sp.all(a == Tinv)

    def test_conduit_conductance_medium(self):
        self.phys.add_model(propname='throat.conduit_conductance',
                            throat_conductance='throat.diffusive_conductance',
                            model=pm.multiphase.conduit_conductance,
                            mode='medium', factor=0)
        Tinv = [1, 2, 3, 4]
        a = sp.where(self.phase['throat.conduit_conductance'] == 0)[0]
        assert sp.all(a == Tinv)

    def test_conduit_conductance_loose(self):
        self.phys.add_model(propname='throat.conduit_conductance',
                            throat_conductance='throat.diffusive_conductance',
                            model=pm.multiphase.conduit_conductance,
                            mode='loose', factor=0)
        Tinv = [1, 2, 3]
        a = sp.where(self.phase['throat.conduit_conductance'] == 0)[0]
        assert sp.all(a == Tinv)

    def test_late_throat_filling(self):
        self.phase['throat.pc_star'] = 1000
        self.phase['throat.pressure'] = 1000
        self.phys.add_model(propname='throat.nwp_saturation',
                            model=pm.multiphase.late_filling,
                            pressure='throat.pressure',
                            Pc_star='throat.pc_star')
        assert sp.all(self.phase['throat.nwp_saturation'] < 1.0)
        assert sp.all(self.phase['throat.nwp_saturation'] > 0.0)

    def test_late_pore_filling(self):
        self.phase['pore.pc_star'] = 1000
        self.phys.add_model(propname='pore.nwp_saturation',
                            model=pm.multiphase.late_filling,
                            Pc_star='pore.pc_star',
                            pressure='pore.pressure')
        assert sp.all(self.phase['pore.nwp_saturation'] < 1.0)
        assert sp.all(self.phase['pore.nwp_saturation'] > 0.0)


if __name__ == '__main__':

    t = MultiPhaseModelsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

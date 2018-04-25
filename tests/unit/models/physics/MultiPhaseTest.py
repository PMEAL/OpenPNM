import openpnm as op
import scipy as sp
from numpy.testing import assert_approx_equal
import openpnm.models.physics as pm


class MultiPhaseTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        Ps = self.net.Ps
        Ts = self.net.Ts
        self.geo = op.geometry.GenericGeometry(network=self.net, pores=Ps,
                                               throats=Ts)
        self.phase = op.phases.GenericPhase(network=self.net)
        self.phase['pore.occupancy'] = sp.ones(self.net.Np)
        self.phase['throat.occupancy'] = sp.ones(self.net.Nt)
        self.phase['pore.occupancy'] = 0.
        self.phase['throat.occupancy'] = 0.
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.capillary_pressure'] = 5000
        self.phys['throat.diffusive_conductance'] = 2

#    def test_conduit_conductance_strict(self):
#        self.phase['pore.occupancy'][[19, 25]] = 0
#        t1 = self.net.Ts[self.phase['throat.occupancy'] == 0]
#        t2 = self.net.Ts[~sp.in1d(self.net.Ts, t1)]
#        self.phys.add_model(propname='throat.cond_conductance',
#                            throat_conductance='throat.diffusive_conductance',
#                            model=pm.multiphase.conduit_conductance,
#                            mode='strict', factor=0)
#        assert ~sp.all(self.phase['throat.cond_conductance'] == 0)
#        assert sp.all(self.phase['throat.cond_conductance'][t1] == 0)
#        assert ~sp.all(self.phase['throat.cond_conductance'][t2] != 0)
#
#    def test_conduit_conductance_medium(self):
#        self.phase['pore.occupancy'][[19, 25]] = 0
#        self.phase['pore.occupancy'][6] = 0
#        t1 = self.net.Ts[self.phase['throat.occupancy'] == 0.]
#        t2 = self.net.find_connecting_throat(6, 7)[0]
#        t3 = self.net.Ts[~sp.in1d(self.net.Ts, sp.concatenate((t1, t2)))]
#        self.phys.add_model(propname='throat.cond_conductance',
#                            throat_conductance='throat.diffusive_conductance',
#                            model=pm.multiphase.conduit_conductance,
#                            mode='medium', factor=0)
#        assert sp.all(self.phase['throat.cond_conductance'][t1] == 0)
#        assert self.phase['throat.cond_conductance'][t2] == 0
#        assert sp.all(self.phase['throat.cond_conductance'][t3] != 0)
#
#    def test_conduit_conductance_loose(self):
#        self.phase['pore.occupancy'][[19, 20]] = 0
#        t1 = self.net.Ts[self.phase['throat.occupancy'] == 0]
#        t2 = self.net.Ts[~sp.in1d(self.net.Ts, t1)]
#        self.phys.add_model(propname='throat.cond_conductance',
#                            throat_conductance='throat.diffusive_conductance',
#                            model=pm.multiphase.conduit_conductance,
#                            mode='loose', factor=0)
#        assert sp.all(self.phase['throat.cond_conductance'][t1] == 0)
#        assert sp.all(self.phase['throat.cond_conductance'][t2] != 0)
#
#    def test_late_throat_filling(self):
#        pass
#
#    def test_late_pore_filling(self):
#        self.phys.add_model(propname='pore.late_p',
#                            model=pm.multiphase.late_pore_filling,
#                            Pc=5500)
#        p = self.phase['pore.occupancy'] > 0
#        assert_approx_equal(self.phase['pore.late_p'][p].mean(),
#                            0.84973704)


if __name__ == '__main__':

    t = MultiPhaseTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

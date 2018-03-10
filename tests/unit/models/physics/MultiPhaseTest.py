import OpenPNM
import scipy as sp
import OpenPNM.Physics.models as pm


class MultiPhaseTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phase['pore.occupancy'] = sp.ones(self.net.Np)
        self.phase['throat.occupancy'] = sp.ones(self.net.Nt)
        p = sp.arange(7, 10)
        t = self.net.find_neighbor_throats(p)
        self.phase['pore.occupancy'][p] = 0.
        self.phase['throat.occupancy'][t] = 0.
        Ps1 = sp.arange(10)
        Ts1 = self.net.find_neighbor_throats(Ps1)
        self.phys1 = OpenPNM.Physics.GenericPhysics(network=self.net,
                                                    phase=self.phase,
                                                    pores=Ps1, throats=Ts1)
        Ps2 = sp.arange(10, self.net.Np)
        Ts2 = self.net.Ts[~sp.in1d(self.net.Ts, Ts1)]
        self.phys2 = OpenPNM.Physics.GenericPhysics(network=self.net,
                                                    phase=self.phase,
                                                    pores=Ps2, throats=Ts2)
        self.phys1['throat.capillary_pressure'] = 5000
        self.phys1['throat.diffusive_conductance'] = 2
        self.phys2['throat.capillary_pressure'] = 5000
        self.phys2['throat.diffusive_conductance'] = 3

    def test_conduit_conductance_strict(self):
        self.phase['pore.occupancy'][[19, 25]] = 0
        t1 = self.net.Ts[self.phase['throat.occupancy'] == 0]
        t2 = self.net.Ts[~sp.in1d(self.net.Ts, t1)]
        self.phys1.models.add(propname='throat.cond_conductance',
                              throat_conductance='throat.diffusive_conductance',
                              model=pm.multiphase.conduit_conductance,
                              mode='strict', factor=0)
        self.phys2.models.add(propname='throat.cond_conductance',
                              throat_conductance='throat.diffusive_conductance',
                              model=pm.multiphase.conduit_conductance,
                              mode='strict', factor=0)
        assert ~sp.all(self.phase['throat.cond_conductance'] == 0)
        assert sp.all(self.phase['throat.cond_conductance'][t1] == 0)
        assert ~sp.all(self.phase['throat.cond_conductance'][t2] != 0)

    def test_conduit_conductance_medium(self):
        self.phase['pore.occupancy'][[19, 25]] = 0
        self.phase['pore.occupancy'][6] = 0
        t1 = self.net.Ts[self.phase['throat.occupancy'] == 0.]
        t2 = self.net.find_connecting_throat(6, 7)[0]
        t3 = self.net.Ts[~sp.in1d(self.net.Ts, sp.concatenate((t1, t2)))]
        self.phys1.models.add(propname='throat.cond_conductance',
                              throat_conductance='throat.diffusive_conductance',
                              model=pm.multiphase.conduit_conductance,
                              mode='medium', factor=0)
        self.phys2.models.add(propname='throat.cond_conductance',
                              throat_conductance='throat.diffusive_conductance',
                              model=pm.multiphase.conduit_conductance,
                              mode='medium', factor=0)
        assert sp.all(self.phase['throat.cond_conductance'][t1] == 0)
        assert self.phase['throat.cond_conductance'][t2] == 0
        assert sp.all(self.phase['throat.cond_conductance'][t3] != 0)

    def test_conduit_conductance_loose(self):
        self.phase['pore.occupancy'][[19, 20]] = 0
        t1 = self.net.Ts[self.phase['throat.occupancy'] == 0]
        t2 = self.net.Ts[~sp.in1d(self.net.Ts, t1)]
        self.phys1.models.add(propname='throat.cond_conductance',
                              throat_conductance='throat.diffusive_conductance',
                              model=pm.multiphase.conduit_conductance,
                              mode='loose', factor=0)
        self.phys2.models.add(propname='throat.cond_conductance',
                              throat_conductance='throat.diffusive_conductance',
                              model=pm.multiphase.conduit_conductance,
                              mode='loose', factor=0)
        assert sp.all(self.phase['throat.cond_conductance'][t1] == 0)
        assert sp.all(self.phase['throat.cond_conductance'][t2] != 0)

    def test_late_throat_filling(self):
        pass

    def test_late_pore_filling(self):
        self.phys1.models.add(propname='pore.late_p',
                              model=pm.multiphase.late_pore_filling,
                              Pc=5500)
        self.phys2.models.add(propname='pore.late_p',
                              model=pm.multiphase.late_pore_filling,
                              Pc=5500)
        p = self.phase['pore.occupancy'] > 0
        assert sp.allclose(self.phase['pore.late_p'][p], 0.84973704)

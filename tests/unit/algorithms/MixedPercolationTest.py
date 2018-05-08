import openpnm as op
import numpy as np
from openpnm.algorithms import MixedInvasionPercolation as mp
import matplotlib.pyplot as plt
import openpnm.models as mods


plt.close('all')
wrk = op.Workspace()
wrk.loglevel = 50


class MixedPercolationTest:

    def setup_class(self):
        # Create Topological Network object
        self.net = op.network.Cubic([5,5,1], spacing=1)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.pores(),
                                               throats=self.net.throats())
        self.geo['pore.diameter'] = 0.5
        self.geo['throat.diameter'] = 0.25
        self.geo.add_model(propname='throat.length',
                           model=mods.geometry.throat_length.straight,
                           L_negative=1e-12,
                           pore_diameter='pore.diameter')
        self.geo.add_model(propname='throat.volume',
                           model=mods.geometry.throat_volume.cylinder,
                           throat_diameter='throat.diameter',
                           throat_length='throat.length')
        self.geo.add_model(propname='pore.volume',
                           model=mods.geometry.pore_volume.sphere,
                           pore_diameter='pore.diameter')
        self.phase = op.phases.Air(network=self.net)
        self.def_phase = op.phases.Water(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.inlets = [0]
        self.outlets = [24]

    def run_mp(self, trapping, partial, snap, plot=False, flowrate=None):
        IP_1 = mp(network=self.net)
        IP_1.settings['partial_saturation']=partial
        IP_1.settings['snap_off']=snap
        IP_1.setup(phase=self.phase,
                   def_phase=self.def_phase)
        IP_1.set_inlets(pores=self.inlets)    
        IP_1.run()
        if trapping:
            IP_1.set_outlets(self.outlets)
            IP_1.apply_trapping()
        inv_points = np.arange(0, 100, 1)
        # returns data as well as plotting
        alg_data = IP_1.plot_drainage_curve(inv_points=inv_points,
                                            lpf=False)
        IP_1.return_results()
        if plot:
            plt.figure()
            plt.imshow(self.phase['pore.invasion_sequence'].reshape([5, 5]),
                       cmap=plt.get_cmap('Blues'))
        else:
            plt.close()
        if flowrate is not None:
            IP_1.apply_flow(flowrate=flowrate)
        return alg_data

    def test_case_throats_sequential(self):
        # Throats only
        # Sequential
        net = self.net
        phys = self.phys
        phys['throat.capillary_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.capillary_pressure']=0.0
        dat_a = self.run_mp(False, False, False)
        # Sequential w. trapping
        dat_b = self.run_mp(True, False, False)
        assert np.all(dat_a[1]==dat_b[1])

    def test_case_throats_random(self):
        # Throats only
        # Random
        net = self.net
        phys = self.phys
        np.random.seed(2)
        phys['throat.capillary_pressure']=np.random.random(net.Nt)*net.Nt
        phys['pore.capillary_pressure']=0.0
        dat_c = self.run_mp(False, False, False)
        # Random w. trapping 
        np.random.seed(2)
        dat_d = self.run_mp(True, False, False)
        assert np.all(dat_d[1]<=dat_c[1])

    def test_case_pores_sequential(self):        
        # Pores only
        # Sequential
        net = self.net
        phys = self.phys
        phys['throat.capillary_pressure']=0.0
        phys['pore.capillary_pressure']=np.arange(0, net.Np, dtype=float)
        dat_e = self.run_mp(False, False, False)
        # Sequential w. trapping
        dat_f = self.run_mp(True, False, False)
        assert np.all(dat_e[1]==dat_f[1])

    def test_case_pores_random(self):
        # Random
        net = self.net
        phys = self.phys
        np.random.seed(2)
        phys['throat.capillary_pressure']=0.0
        phys['pore.capillary_pressure']=np.random.random(net.Np)*net.Np
        dat_g = self.run_mp(False, False, False)
        # Random w. trapping 
        np.random.seed(2)
        dat_h = self.run_mp(True, False, False)
        assert np.all(dat_h[1]<=dat_g[1])

    def test_case_mixed_sequential(self): 
        # Pores and Throats
        # Sequential
        net = self.net
        phys = self.phys
        phys['throat.capillary_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.capillary_pressure']=np.arange(0, net.Np, dtype=float)
        dat_i = self.run_mp(False, False, False)
        # Sequential w. trapping
        dat_j = self.run_mp(True, False, False)
        assert np.all(dat_i[1]==dat_j[1])

    def test_case_mixed_random(self):        
        # Random
        net = self.net
        phys = self.phys
        np.random.seed(2)
        phys['throat.capillary_pressure']=np.random.random(net.Nt)*net.Nt
        phys['pore.capillary_pressure']=np.random.random(net.Np)*net.Np
        dat_k = self.run_mp(False, False, False)
        # Random w. trapping 
        np.random.seed(2)
        dat_l = self.run_mp(True, False, False)
        assert np.all(dat_l[1]<=dat_k[1])


    def test_snap_off(self):
        # Throats only
        # Sequential
        net = self.net
        phys = self.phys
        phys['throat.capillary_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.capillary_pressure']=0.0
        dat_m = self.run_mp(False, False, False)
        # Sequential w. snap-off
        phys['throat.snap_off']=100.0 # This pressure is higher than burst
        T = 10
        [P1, P2] = self.net['throat.conns'][T]
        phys['throat.snap_off'][T]=0.5 # This pressure is lower than burst
        dat_n = self.run_mp(False, False, True)
        assert self.phase['pore.invasion_pressure'][P1] == 0.5
        assert self.phase['pore.invasion_pressure'][P2] == 0.5
        assert self.phase['throat.invasion_pressure'][T] == 0.5
        assert ~np.all(dat_m[1]-dat_n[1]==0)

    def test_partial(self):
        # Throats only
        # Sequential
        net = self.net
        phys = self.phys
        phys['throat.capillary_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.capillary_pressure']=0.0
        dat_o = self.run_mp(False, False, False)
        # Sequential w. partial
        T = 10
        [P1, P2] = self.net['throat.conns'][T]
        self.phase['pore.occupancy'] = False
        self.phase['throat.occupancy'] = False
        self.phase['pore.occupancy'][P1] = True
        self.phase['pore.occupancy'][P2] = True
        dat_p = self.run_mp(False, True, False, True)
        assert self.phase['pore.invasion_pressure'][P1] == -np.inf
        assert self.phase['pore.invasion_pressure'][P2] == -np.inf
        assert self.phase['throat.invasion_pressure'][T] == -np.inf
        assert ~np.all(dat_o[1]-dat_p[1]==0)

    def test_apply_flow_rate(self):
        t = self
        pvol = np.sum(t.net['pore.volume'])
        tvol = np.sum(t.net['throat.volume'])
        tot = pvol+tvol
        net = self.net
        phys = self.phys
        phys['throat.capillary_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.capillary_pressure']=0.0
        self.run_mp(False, False, False, flowrate=tot)
        assert 'throat.invasion_time' in self.phase.props()

if __name__ == '__main__':
    t = MixedPercolationTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

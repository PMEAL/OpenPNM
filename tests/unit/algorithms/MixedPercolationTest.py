import openpnm as op
import numpy as np
from openpnm.algorithms import MixedInvasionPercolation as mp
import matplotlib.pyplot as plt
import openpnm.models as mods


plt.close('all')
wrk = op.Workspace()
wrk.loglevel = 50


class MixedPercolationTest:

    def setup_class(self, Np=5):
        # Create Topological Network object
        self.net = op.network.Cubic([Np, Np, 1], spacing=1)
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
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.inlets = [0]
        self.outlets = [Np*Np - 1]

    def run_mp(self, trapping=False, residual=False, snap=False,
               plot=False, flowrate=None):
        IP_1 = mp(network=self.net)
        IP_1.settings['snap_off']=snap
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=self.inlets)
        if residual:
            IP_1.set_residual(pores=self.phase['pore.occupancy'])
        IP_1.run()
        if trapping:
            IP_1.set_outlets(self.outlets)
            IP_1.apply_trapping()
        inv_points = np.arange(0, 100, 1)
        # returns data as well as plotting
        alg_data = IP_1.plot_drainage_curve(inv_points=inv_points,
                                            lpf=False)
        IP_1.results()
        if plot:
            plt.figure()
            l = np.sqrt(self.net.Np).astype(int)
            plt.imshow(IP_1['pore.invasion_sequence'].reshape([l, l]),
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
        phys['throat.entry_pressure'] = np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure'] = 0.0
        dat_a = self.run_mp(False, False, False)
        # Sequential w. trapping
        dat_b = self.run_mp(True, False, False)
#        assert np.all(dat_a[1] == dat_b[1])

    def test_case_throats_random(self):
        # Throats only
        # Random
        net = self.net
        phys = self.phys
        np.random.seed(2)
        phys['throat.entry_pressure'] = np.random.random(net.Nt)*net.Nt
        phys['pore.entry_pressure'] = 0.0
        dat_c = self.run_mp(False, False, False)
        # Random w. trapping
        np.random.seed(2)
        dat_d = self.run_mp(True, False, False)
        assert np.all(dat_d[1] <= dat_c[1])

    def test_case_pores_sequential(self):
        # Pores only
        # Sequential
        net = self.net
        phys = self.phys
        phys['throat.entry_pressure'] = 0.0
        phys['pore.entry_pressure'] = np.arange(0, net.Np, dtype=float)
        dat_e = self.run_mp(False, False, False)
        # Sequential w. trapping
        dat_f = self.run_mp(True, False, False)
#        assert np.all(dat_e[1] == dat_f[1])

    def test_case_pores_random(self):
        # Random
        net = self.net
        phys = self.phys
        np.random.seed(2)
        phys['throat.entry_pressure'] = 0.0
        phys['pore.entry_pressure'] = np.random.random(net.Np)*net.Np
        dat_g = self.run_mp(False, False, False)
        # Random w. trapping
        np.random.seed(2)
        dat_h = self.run_mp(True, False, False)
        assert np.all(dat_h[1] <= dat_g[1])

    def test_case_mixed_sequential(self):
        # Pores and Throats
        # Sequential
        net = self.net
        phys = self.phys
        phys['throat.entry_pressure'] = np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure'] = np.arange(0, net.Np, dtype=float)
        dat_i = self.run_mp(False, False, False)
        # Sequential w. trapping
        dat_j = self.run_mp(True, False, False)
#        assert np.all(dat_i[1] == dat_j[1])

    def test_case_mixed_random(self):
        # Random
        net = self.net
        phys = self.phys
        np.random.seed(2)
        phys['throat.entry_pressure']=np.random.random(net.Nt)*net.Nt
        phys['pore.entry_pressure']=np.random.random(net.Np)*net.Np
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
        phys['throat.entry_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure']=0.0
        dat_m = self.run_mp(False, False, False)
        # Sequential w. snap-off
        phys['throat.snap_off']=100.0  # This pressure is higher than burst
        T = 10
        [P1, P2] = self.net['throat.conns'][T]
        phys['throat.snap_off'][T]=0.5  # This pressure is lower than burst
        dat_n = self.run_mp(False, False, True)
        assert self.phase['pore.invasion_pressure'][P1] == 0.5
        assert self.phase['pore.invasion_pressure'][P2] == 0.5
        assert self.phase['throat.invasion_pressure'][T] == 0.5
        assert ~np.all(dat_m[1]-dat_n[1]==0)

    def test_residual(self):
        # Throats only
        # Sequential
        net = self.net
        phys = self.phys
        phys['throat.entry_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure']=0.0
        dat_o = self.run_mp(False, False, False)
        # Sequential w. partial
        T = 10
        [P1, P2] = self.net['throat.conns'][T]
        self.phase['pore.occupancy'] = False
        self.phase['throat.occupancy'] = False
        self.phase['pore.occupancy'][P1] = True
        self.phase['pore.occupancy'][P2] = True
        dat_p = self.run_mp(False, True, False, False)
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
        phys['throat.entry_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure']=0.0
        self.run_mp(False, False, False, flowrate=tot)
        assert 'throat.invasion_time' in self.phase.props()

    def test_max_pressure(self):
        net = self.net
        phys = self.phys
        phys['throat.entry_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure']=0.0
        IP_1 = mp(network=self.net)
        IP_1.settings['partial_saturation']=False
        IP_1.settings['snap_off']=False
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=self.inlets)
        IP_1.run(max_pressure=20)
        IP_1.results()
        inv_Pc = self.phase['pore.invasion_pressure']
        inv_Pc = inv_Pc[~np.isinf(inv_Pc)]
        assert inv_Pc.max() <= 20

    def test_drainage_curve(self):
        net = self.net
        phys = self.phys
        phys['throat.entry_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure']=0.0
        IP_1 = mp(network=self.net)
        self.phase['pore.occupancy'] = False
        self.phase['throat.occupancy'] = False
        IP_1.settings['snap_off']=False
        IP_1.setup(phase=self.phase)
        inv_points = np.arange(0, 100, 1, dtype=float)
        sat = np.zeros_like(inv_points)
        tot_vol = (np.sum(self.net['pore.volume']) +
                   np.sum(self.net['throat.volume']))
        for i, Pc in enumerate(inv_points):
            IP_1.reset()
            IP_1.set_inlets(pores=self.inlets)
            IP_1.set_residual(pores=self.phase['pore.occupancy'])
            IP_1.run(max_pressure=Pc)
            IP_1.results()
            Pinv_Pc = self.phase['pore.invasion_pressure']
            Tinv_Pc = self.phase['throat.invasion_pressure']
            sat[i] += np.sum(self.net['pore.volume'][Pinv_Pc<np.inf])
            sat[i] += np.sum(self.net['throat.volume'][Tinv_Pc<np.inf])
        assert sat.max()/tot_vol == 1.0

    def test_cluster_merging(self):
        phys = self.phys
        phys['throat.entry_pressure']=0.0
        Pc = np.array([[0.0, 1.0, 2.0, 1.0, 0.0],
                       [3.0, 4.0, 5.0, 4.0, 3.0],
                       [6.0, 7.0, 8.0, 7.0, 6.0],
                       [9.0, 10.0, 11.0, 10.0, 9.0],
                       [12.0, 13.0, 14.0, 13.0, 12.0]])
        phys['pore.entry_pressure']=Pc.flatten()

        IP_1 = mp(network=self.net)
        IP_1.settings['partial_saturation']=False
        IP_1.settings['snap_off']=False
        IP_1.setup(phase=self.phase)
        # Set the inlets as the pores with zero entry Pc
        IP_1.set_inlets(clusters=[[0], [4]])
        IP_1.run()
        IP_1.results()
        # Clusters should merge on first row and all pores after the first row
        # should be part of the same cluster
        assert len(np.unique(self.phase['pore.cluster'][5:])) == 1

    def test_connected_residual_clusters(self):
        net = self.net
        phys = self.phys
        phys['throat.entry_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure']=np.arange(0, net.Np, dtype=float)
        IP_1 = mp(network=self.net)
        IP_1.settings['residual_saturation']=True
        IP_1.settings['snap_off']=False
        IP_1.setup(phase=self.phase)
        T = 20
        [P1, P2] = self.net['throat.conns'][T]
        self.phase['pore.occupancy'] = False
        self.phase['throat.occupancy'] = False
        self.phase['pore.occupancy'][P1] = True
        self.phase['pore.occupancy'][P2] = True
        IP_1.set_inlets(pores=self.inlets)
        assert len(IP_1.queue) == 1

    def test_disconnected_residual_clusters(self):
        net = self.net
        phys = self.phys
        phys['throat.entry_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure']=np.arange(0, net.Np, dtype=float)
        IP_1 = mp(network=self.net)
        IP_1.settings['snap_off']=False
        IP_1.setup(phase=self.phase)
        T = 20
        [P1, P2] = self.net['throat.conns'][T]
        self.phase['pore.occupancy'] = False
        self.phase['throat.occupancy'] = False
        self.phase['pore.occupancy'][P1] = False
        self.phase['pore.occupancy'][P2] = True
        IP_1.set_inlets(pores=self.inlets)
        IP_1.set_residual(pores=self.phase['pore.occupancy'])
        assert len(IP_1.queue) == 2

    def test_big_clusters(self):
        self.setup_class(Np=10)
        net = self.net
        phys = self.phys
        np.random.seed(1)
        phys['throat.entry_pressure'] = 0.0
        phys['pore.entry_pressure'] = np.random.random(net.Np)*net.Np
        self.inlets = net.pores('left')
        self.outlets = None
        np.random.seed(1)
        self.phase['pore.occupancy'] = False
        self.phase['throat.occupancy'] = False
        self.phase['pore.occupancy'] = np.random.random(net.Np) < 0.25
        IP_1 = mp(network=self.net)
        IP_1.settings['snap_off'] = False
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=self.inlets)
        IP_1.set_residual(pores=self.phase['pore.occupancy'])
        IP_1.run()
        IP_1.results()
        assert np.all(self.phase['pore.invasion_sequence'] > -1)
        assert len(np.unique(self.phase['pore.cluster'])) > 1

    def test_big_clusters_trapping(self):
        self.setup_class(Np=10)
        net = self.net
        phys = self.phys
        np.random.seed(1)
        phys['throat.entry_pressure'] = 0.0
        phys['pore.entry_pressure'] = np.random.random(net.Np)*net.Np
        self.inlets = net.pores('left')
        self.outlets = net.pores('right')
        np.random.seed(1)
        self.phase['pore.occupancy'] = False
        self.phase['throat.occupancy'] = False
        self.phase['pore.occupancy'] = np.random.random(net.Np) < 0.25
        IP_1 = mp(network=self.net)
        IP_1.settings['snap_off'] = False
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=self.inlets)
        IP_1.set_residual(pores=self.phase['pore.occupancy'])
        IP_1.run()
        IP_1.set_outlets(self.outlets)
        IP_1.apply_trapping()
        IP_1.results()
#        assert np.sum(IP_1['pore.trapped']) ==35

    def test_invade_isolated_Ts(self):
        self.setup_class(Np=10)
        net = self.net
        phys = self.phys
        np.random.seed(1)
        phys['throat.entry_pressure'] = 0.0
        phys['pore.entry_pressure'] = np.random.random(net.Np)*net.Np
        self.inlets = net.pores('left')
        self.outlets = None
        IP_1 = mp(network=self.net)
        IP_1.settings['invade_isolated_Ts'] = False
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=self.inlets)
        IP_1.run()
        IP_1.results()
        save_seq = IP_1['throat.invasion_sequence'].copy()
        IP_1.settings['invade_isolated_Ts'] = True
        IP_1.reset()
        IP_1.set_inlets(pores=self.inlets)
        IP_1.run()
        IP_1.results()
        assert np.any(IP_1['throat.invasion_sequence']-save_seq != 0)

    def test_terminate_clusters(self):
        self.setup_class(Np=10)
        net = self.net
        phys = self.phys
        np.random.seed(1)
        phys['throat.entry_pressure'] = 0.0
        phys['pore.entry_pressure'] = np.random.random(net.Np)*net.Np
        inlets = net.pores('left')
        outlets = net.pores('right')
        IP_1 = mp(network=self.net)
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=inlets)
        IP_1.set_outlets(pores=outlets)
        IP_1.run()
        IP_1.results()
        assert np.any(IP_1['throat.invasion_sequence'][outlets] > -1)
        assert np.any(IP_1['throat.invasion_sequence'] == -1)


if __name__ == '__main__':
    t = MixedPercolationTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

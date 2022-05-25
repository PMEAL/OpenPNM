import numpy as np
import openpnm as op
import openpnm.models.geometry as gm
from openpnm.algorithms import MixedInvasionPercolation as mp
import matplotlib.pyplot as plt


plt.close('all')
ws = op.Workspace()
ws.loglevel = 50


class MixedPercolationTest:

    def setup_class(self, Np=5):
        ws.clear()
        # Create Topological Network object
        self.net = op.network.Cubic([Np, Np, 1], spacing=1)
        self.net.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders()
        )
        self.net.regenerate_models()
        self.net['pore.diameter'] = 0.5
        self.net['throat.diameter'] = 0.25
        self.phase = op.phase.Air(network=self.net)
        self.inlets = [0]
        self.outlets = [Np*Np - 1]

    def run_mp(self, trapping=False, residual=False, snap=False,
               plot=False, flowrate=None):
        IP_1 = mp(network=self.net)
        if snap:
            IP_1.settings['snap_off'] = 'throat.snap_off'
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
        alg_data = IP_1.get_intrusion_data(inv_points=inv_points)
        self.phase.update(IP_1.results(Pc=inv_points.max()))
        if plot:
            plt.figure()
            L = np.sqrt(self.net.Np).astype(int)
            plt.imshow(IP_1['pore.invasion_sequence'].reshape([L, L]),
                       cmap=plt.get_cmap('Blues'))
        if flowrate is not None:
            IP_1.apply_flow(flowrate=flowrate)
        self.alg = IP_1
        return alg_data

    def test_case_throats_sequential(self):
        # Throats only
        # Sequential
        self.phase['throat.entry_pressure'] = np.arange(0, self.net.Nt, dtype=float)
        self.phase['pore.entry_pressure'] = 0.0
        dat_a = self.run_mp(False, False, False)
        # Sequential w. trapping
        dat_b = self.run_mp(True, False, False)
        assert np.all(dat_a.S_tot == dat_b.S_tot)

    def test_case_throats_random(self):
        # Throats only
        # Random
        np.random.seed(2)
        self.phase['throat.entry_pressure'] = \
            np.random.random(self.net.Nt)*self.net.Nt
        self.phase['pore.entry_pressure'] = 0.0
        dat_c = self.run_mp(False, False, False)
        # Random w. trapping
        np.random.seed(2)
        dat_d = self.run_mp(True, False, False)
        assert np.all(dat_d.S_tot <= dat_c.S_tot)

    def test_case_pores_sequential(self):
        # Pores only
        # Sequential
        self.phase['throat.entry_pressure'] = 0.0
        self.phase['pore.entry_pressure'] = np.arange(0, self.net.Np, dtype=float)
        dat_e = self.run_mp(False, False, False)
        # Sequential w. trapping
        dat_f = self.run_mp(True, False, False)
        assert np.all(dat_e.S_tot == dat_f.S_tot)

    def test_case_pores_random(self):
        # Random
        np.random.seed(2)
        self.phase['throat.entry_pressure'] = 0.0
        self.phase['pore.entry_pressure'] = np.random.random(self.net.Np)*self.net.Np
        dat_g = self.run_mp(False, False, False)
        # Random w. trapping
        np.random.seed(2)
        dat_h = self.run_mp(True, False, False)
        assert np.all(dat_h.S_tot <= dat_g.S_tot)

    def test_case_mixed_sequential(self):
        # Pores and Throats
        # Sequential
        self.phase['throat.entry_pressure'] = np.arange(0, self.net.Nt, dtype=float)
        self.phase['pore.entry_pressure'] = np.arange(0, self.net.Np, dtype=float)
        dat_i = self.run_mp(False, False, False)
        # Sequential w. trapping
        dat_j = self.run_mp(True, False, False)
        assert np.all(dat_i.S_tot==dat_j.S_tot)

    def test_case_mixed_random(self):
        # Random
        np.random.seed(2)
        self.phase['throat.entry_pressure'] = \
            np.random.random(self.net.Nt)*self.net.Nt
        self.phase['pore.entry_pressure'] = \
            np.random.random(self.net.Np)*self.net.Np
        dat_k = self.run_mp(False, False, False)
        # Random w. trapping
        np.random.seed(2)
        dat_l = self.run_mp(True, False, False)
        assert np.all(dat_l.S_tot <= dat_k.S_tot)

    def test_snap_off(self):
        # Throats only
        # Sequential
        net = self.net
        phys = self.phase
        phys['throat.entry_pressure'] = np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure'] = 0.0
        dat_m = self.run_mp(False, False, False)
        # Sequential w. snap-off
        phys['throat.snap_off'] = 100.0  # This pressure is higher than burst
        T = 10
        [P1, P2] = self.net['throat.conns'][T]
        phys['throat.snap_off'][T] = 0.5  # This pressure is lower than burst
        dat_n = self.run_mp(False, False, True)
        assert self.alg['pore.invasion_pressure'][P1] == 0.5
        assert self.alg['pore.invasion_pressure'][P2] == 0.5
        assert self.alg['throat.invasion_pressure'][T] == 0.5
        assert ~np.all(dat_m.S_tot-dat_n.S_tot == 0)

    def test_residual(self):
        # Throats only
        # Sequential
        net = self.net
        phys = self.phase
        phys['throat.entry_pressure'] = np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure'] = 0.0
        dat_o = self.run_mp(False, False, False)
        # Sequential w. partial
        T = 10
        [P1, P2] = self.net['throat.conns'][T]
        self.phase['pore.occupancy'] = False
        self.phase['throat.occupancy'] = False
        self.phase['pore.occupancy'][P1] = True
        self.phase['pore.occupancy'][P2] = True
        dat_p = self.run_mp(False, True, False, False)
        assert self.alg['pore.invasion_pressure'][P1] == -np.inf
        assert self.alg['pore.invasion_pressure'][P2] == -np.inf
        assert self.alg['throat.invasion_pressure'][T] == -np.inf
        assert ~np.all(dat_o.S_tot-dat_p.S_tot == 0)

    def test_apply_flow_rate(self):
        t = self
        pvol = np.sum(t.net['pore.volume'])
        tvol = np.sum(t.net['throat.volume'])
        tot = pvol+tvol
        net = self.net
        phys = self.phase
        phys['throat.entry_pressure'] = np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure'] = 0.0
        self.run_mp(False, False, False, flowrate=tot)
        assert 'throat.invasion_time' in self.phase.props()

    def test_max_pressure(self):
        net = self.net
        phys = self.phase
        phys['throat.entry_pressure'] = np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure'] = 0.0
        IP_1 = mp(network=self.net)
        IP_1.settings['partial_saturation'] = False
        IP_1.settings['snap_off'] = ""
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=self.inlets)
        IP_1.run(max_pressure=20)
        IP_1.results(Pc=20)
        inv_Pc = IP_1['pore.invasion_pressure']
        inv_Pc = inv_Pc[~np.isinf(inv_Pc)]
        assert inv_Pc.max() <= 20

    def test_drainage_curve(self):
        net = self.net
        phys = self.phase
        phys['throat.entry_pressure'] = np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure'] = 0.0
        IP_1 = mp(network=self.net)
        self.phase['pore.occupancy'] = False
        self.phase['throat.occupancy'] = False
        IP_1.settings['snap_off'] = ""
        IP_1.setup(phase=self.phase)
        inv_points = np.arange(0, 100, 1, dtype=float)
        sat = np.zeros_like(inv_points)
        tot_vol = self.net['pore.volume'].sum() + self.net['throat.volume'].sum()
        for i, Pc in enumerate(inv_points):
            IP_1.reset()
            IP_1.set_inlets(pores=self.inlets)
            IP_1.set_residual(pores=self.phase['pore.occupancy'])
            IP_1.run(max_pressure=Pc)
            IP_1.results(Pc)
            Pinv_Pc = IP_1['pore.invasion_pressure']
            Tinv_Pc = IP_1['throat.invasion_pressure']
            sat[i] += np.sum(self.net['pore.volume'][Pinv_Pc<np.inf])
            sat[i] += np.sum(self.net['throat.volume'][Tinv_Pc<np.inf])
        assert sat.max()/tot_vol == 1.0

    def test_plot_intrusion_curve(self):
        net = self.net
        phys = self.phase
        phys['throat.entry_pressure']=np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure']=0.0
        self.run_mp(False, False, False)
        fig, ax = plt.subplots()
        self.alg.plot_intrusion_curve(ax=ax)
        plt.close()
        self.alg.plot_intrusion_curve()
        plt.close()

    def test_cluster_merging(self):
        phys = self.phase
        phys['throat.entry_pressure'] = 0.0
        Pc = np.array([[0.0, 1.0, 2.0, 1.0, 0.0],
                       [3.0, 4.0, 5.0, 4.0, 3.0],
                       [6.0, 7.0, 8.0, 7.0, 6.0],
                       [9.0, 10.0, 11.0, 10.0, 9.0],
                       [12.0, 13.0, 14.0, 13.0, 12.0]])
        phys['pore.entry_pressure'] = Pc.flatten()

        IP_1 = mp(network=self.net)
        IP_1.settings['partial_saturation'] = ""
        IP_1.settings['snap_off'] = ""
        IP_1.setup(phase=self.phase)
        # Set the inlets as the pores with zero entry Pc
        IP_1.set_inlets(clusters=[[0], [4]])
        IP_1.run()
        # Clusters should merge on first row and all pores after the first row
        # should be part of the same cluster
        assert len(np.unique(IP_1['pore.cluster'][5:])) == 1

    def test_connected_residual_clusters(self):
        net = self.net
        phys = self.phase
        phys['throat.entry_pressure'] = np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure'] = np.arange(0, net.Np, dtype=float)
        IP_1 = mp(network=self.net)
        IP_1.settings['residual_saturation'] = True
        IP_1.settings['snap_off'] = ""
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
        phys = self.phase
        phys['throat.entry_pressure'] = np.arange(0, net.Nt, dtype=float)
        phys['pore.entry_pressure'] = np.arange(0, net.Np, dtype=float)
        IP_1 = mp(network=self.net)
        IP_1.settings['snap_off'] = ""
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
        phys = self.phase
        np.random.seed(1)
        phys['throat.entry_pressure']=0.0
        phys['pore.entry_pressure']=np.random.random(net.Np)*net.Np
        self.inlets = net.pores('left')
        self.outlets = None
        np.random.seed(1)
        self.phase['pore.occupancy'] = False
        self.phase['throat.occupancy'] = False
        self.phase['pore.occupancy'] = np.random.random(net.Np) < 0.25
        IP_1 = mp(network=self.net)
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=self.inlets)
        IP_1.set_residual(pores=self.phase['pore.occupancy'])
        IP_1.run()
        assert np.all(IP_1['pore.invasion_sequence'] > -1)
        assert len(np.unique(IP_1['pore.cluster'])) > 1

    def test_big_clusters_trapping(self):
        self.setup_class(Np=10)
        net = self.net
        phys = self.phase
        np.random.seed(1)
        phys['throat.entry_pressure']=0.0
        phys['pore.entry_pressure']=np.random.random(net.Np)*net.Np
        self.inlets = net.pores('front')
        self.outlets = net.pores('back')
        np.random.seed(1)
        self.phase['pore.occupancy'] = False
        self.phase['throat.occupancy'] = False
        self.phase['pore.occupancy'] = np.random.random(net.Np) < 0.25
        IP_1 = mp(network=self.net)
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=self.inlets)
        IP_1.set_residual(pores=self.phase['pore.occupancy'])
        IP_1.run()
        IP_1.set_outlets(self.outlets)
        IP_1.apply_trapping()
        assert np.sum(IP_1['pore.trapped'])==35

    def test_invade_isolated_Ts(self):
        self.setup_class(Np=10)
        net = self.net
        phys = self.phase
        np.random.seed(1)
        phys['throat.entry_pressure']=0.0
        phys['pore.entry_pressure']=np.random.random(net.Np)*net.Np
        self.inlets = net.pores('left')
        self.outlets = None
        IP_1 = mp(network=self.net)
        IP_1.settings['invade_isolated_Ts'] = False
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=self.inlets)
        IP_1.run()
        save_seq = IP_1['throat.invasion_sequence'].copy()
        IP_1.settings['invade_isolated_Ts'] = True
        IP_1.reset()
        IP_1.set_inlets(pores=self.inlets)
        IP_1.run()
        assert np.any(IP_1['throat.invasion_sequence']-save_seq != 0)

    def test_terminate_clusters(self):
        self.setup_class(Np=10)
        net = self.net
        phys = self.phase
        np.random.seed(1)
        phys['throat.entry_pressure'] = 0.0
        phys['pore.entry_pressure']=np.random.random(net.Np)*net.Np
        inlets = net.pores('left')
        outlets = net.pores('right')
        IP_1 = mp(network=self.net)
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=inlets)
        IP_1.set_outlets(pores=outlets)
        IP_1.run()
        assert np.any(IP_1['throat.invasion_sequence'][outlets] > -1)
        assert np.any(IP_1['throat.invasion_sequence'] == -1)

    def test_late_filling(self):
        self.setup_class(Np=10)
        net = self.net
        phys = self.phase
        np.random.seed(1)
        phys['throat.entry_pressure'] = np.random.random(net.Nt)*10000 + 5000
        phys['pore.entry_pressure'] = 0.0
        phys.add_model(propname='pore.pc_star',
                       model=op.models.misc.from_neighbor_throats,
                       prop='throat.entry_pressure',
                       mode='min')
        phys.add_model(propname='pore.late_filling',
                       model=op.models.physics.multiphase.late_filling,
                       pressure='pore.pressure',
                       Pc_star='pore.pc_star',
                       eta=1, Swp_star=0.4)
        phys['throat.pc_star'] = phys['throat.entry_pressure']
        phys.add_model(propname='throat.late_filling',
                       model=op.models.physics.multiphase.late_filling,
                       pressure='throat.pressure',
                       Pc_star='throat.pc_star',
                       eta=1, Swp_star=0.2)
        inlets = net.pores('left')
        outlets = net.pores('right')
        IP_1 = mp(network=self.net)
        IP_1.setup(phase=self.phase)
        IP_1.set_inlets(pores=inlets)
        IP_1.set_outlets(pores=outlets)
        IP_1.run()
        inv_points = np.arange(phys['throat.entry_pressure'].min(),
                               phys['throat.entry_pressure'].max(), 100)
        alg_data = IP_1.get_intrusion_data(inv_points=inv_points)
        IP_1.settings['late_pore_filling'] = 'pore.late_filling'
        IP_1.settings['late_throat_filling'] = 'throat.late_filling'
        alg_data_lpf = IP_1.get_intrusion_data(inv_points=inv_points)
        assert np.any(alg_data_lpf.S_tot - alg_data.S_tot < 0.0)
        assert ~np.any(alg_data_lpf.S_tot - alg_data.S_tot > 0.0)

    # def test_bidirectional_entry_pressure(self):
    #     pn = op.network.Cubic(shape=[3, 3, 3], spacing=2.5e-5)
    #     pn.add_model_collection(
    #         op.models.collections.geometry.spheres_and_cylinders()
    #     )
    #     pn.regenerate_models()
    #     pn['throat.diameter'] = 2.0e-5
    #     pn['pore.diameter'] = (np.random.random(pn.Np)+0.5)*1e-5
    #     pn['pore.volume'] = (4/3)*np.pi*(pn['pore.diameter']/2)**3
    #     pn['throat.volume'] = 0.0
    #     pn.add_model(propname='throat.centroid',
    #                  model=op.models.geometry.throat_centroid.pore_coords)
    #     pn.add_model(propname='throat.normal',
    #                  model=op.models.geometry.throat_vector.pore_to_pore)
    #     water = op.phase.Water(network=pn)
    #     water['pore.contact_angle'] = 100
    #     r_tor = 5e-6
    #     pmod = op.models.physics.capillary_pressure.purcell_bidirectional
    #     water.add_model(propname='throat.entry_pressure',
    #                     model=pmod,
    #                     r_toroid=r_tor)
    #     water.add_model(propname='throat.max_pressure',
    #                     model=op.models.physics.meniscus.purcell,
    #                     r_toroid=r_tor,
    #                     mode='max')
    #     water['pore.entry_pressure'] = 0.0
    #     ip = op.algorithms.MixedInvasionPercolation(network=pn)
    #     ip.setup(phase=water)
    #     ip.set_inlets(pores=pn.pores('bottom'))
    #     ip.run()
    #     alg_data = ip.get_intrusion_data()
    #     # Max pressure is all the same but bi-directional touch pressure isn't
    #     # So there will be different invasion points. Using max results in a
    #     # Single invasion point
    #     assert np.any(alg_data.S_pore < 1.0)


if __name__ == '__main__':
    t = MixedPercolationTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

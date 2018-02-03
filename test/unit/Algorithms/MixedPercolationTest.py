import OpenPNM as op
import numpy as np
from OpenPNM.Physics import models as pm
import matplotlib.pyplot as plt
plt.close('all')
wrk = op.Base.Workspace()
wrk.loglevel = 50


class MixedPercolationTest:

    def setup_class(self):
        # Create Topological Network object

        self.fiber_rad = 5e-6
        self.net = op.Network.Delaunay(num_pores=25,
                                       domain_size=[50e-6, 50e-6, 50e-6])
        self.net.add_boundaries()

        Ps = self.net.pores()
        bulk_pores = self.net.pores('boundary', mode='not')
        self.net['pore.bulk'] = False
        self.net['pore.bulk'][bulk_pores] = True
        BTs = self.net.find_neighbor_throats(self.net.pores('boundary'))
        self.net['throat.bulk'] = True
        self.net['throat.bulk'][BTs] = False
        Ts = self.net.find_neighbor_throats(pores=Ps,
                                            mode='intersection',
                                            flatten=True)
        self.geom = op.Geometry.Voronoi(network=self.net,
                                        pores=Ps,
                                        throats=Ts,
                                        fibre_rad=self.fiber_rad,
                                        voxel_vol=True,
                                        vox_len=1e-6)
        self.geom.add_model(propname='throat.volume',
                            model=op.Geometry.models.throat_volume.extrusion)
        pn = self.net
        self.net.trim(throats=pn.throats()[pn["throat.area"] == 0])
        self.net.trim(throats=pn.throats()[pn["throat.perimeter"] == 0])
        self.air = op.Phases.Air(network=pn)
        self.water = op.Phases.Water(network=pn)
        self.water["pore.contact_angle"] = 110
        self.Swp_star = 0.25  # late pore filling
        self.air["pore.contact_angle"] = 70
        self.air["pore.surface_tension"] = self.water["pore.surface_tension"]
        self.inv_points = np.linspace(0, 20000, 21)

    def process_physics(self, model='purcell', snap_off=True):
        # Clean up already made phys
        try:
            phys_air = self.net.physics('mp_phys_a')[0]
            wrk.purge_object(phys_air)
        except:
            pass
        try:
            phys_water = self.net.physics('mp_phys_w')[0]
            wrk.purge_object(phys_water)
        except:
            pass
        Ps = self.net.pores()
        Ts = self.net.throats()
        phys_air = op.Physics.Standard(network=self.net,
                                       phase=self.air,
                                       pores=Ps,
                                       throats=Ts,
                                       name='mp_phys_a')
        phys_water = op.Physics.Standard(network=self.net,
                                         phase=self.water,
                                         pores=Ps,
                                         throats=Ts,
                                         name='mp_phys_w')
        throat_diam = 'throat.diameter'
        pore_diam = 'pore.indiameter'
        if model == 'purcell':
            pmod = pm.capillary_pressure.purcell
        elif model == 'purcell_bi':
            pmod = pm.capillary_pressure.purcell_bi
        elif model == 'sinusoidal':
            pass
        phys_water.models.add(propname='throat.capillary_pressure',
                              model=pmod,
                              r_toroid=self.fiber_rad,
                              diameter=throat_diam,
                              h_max=pore_diam)
        phys_air.models.add(propname='throat.capillary_pressure',
                            model=pmod,
                            r_toroid=self.fiber_rad,
                            diameter=throat_diam,
                            h_max=pore_diam)
        if snap_off:
            phys_air.models.add(propname='throat.snap_off',
                                model=pm.capillary_pressure.ransohoff_snap_off,
                                throat_diameter=throat_diam,
                                wavelength=self.fiber_rad)
            phys_air['throat.snap_off'] = np.abs(phys_air['throat.snap_off'])
        phys_air['pore.capillary_pressure'] = 0
        phys_water['pore.capillary_pressure'] = 0
        BPs = self.net.pores('boundary')
        NBPs = self.net.find_neighbor_pores(BPs, flatten=False)
        boundary_neighbors = []
        for NBP in NBPs:
            boundary_neighbors.append(NBP[0])
        NBPs = np.asarray(boundary_neighbors)
        wPc_NBPs = phys_water["pore.capillary_pressure"][NBPs]
        phys_water["pore.capillary_pressure"][BPs] = wPc_NBPs
        aPc_NBPs = phys_air["pore.capillary_pressure"][NBPs]
        phys_air["pore.capillary_pressure"][BPs] = aPc_NBPs

    def run_alg(self, inv_phase=None, def_phase=None,
                lpf=False,
                snap_off=False,
                coop_fill=False,
                partial=False,
                trapping=False,
                cap_model='purcell'):
        self.process_physics(model=cap_model, snap_off=snap_off)
        inlets = self.net.pores(labels=['bottom_boundary'])
        outlets = self.net.pores(labels=['top_boundary'])
        in_step = 2
        ip_inlets = [inlets[x] for x in range(0, len(inlets), in_step)]
        inlet_inv_seq = -1

        IP_1 = op.Algorithms.MixedPercolation(network=self.net)
        IP_1.setup(phase=inv_phase,
                   def_phase=def_phase,
                   inlets=ip_inlets,
                   inlet_inv_seq=inlet_inv_seq,
                   snap_off=snap_off,
                   partial=partial)
        if coop_fill:
            IP_1.setup_coop_filling(capillary_model=cap_model,
                                    inv_points=self.inv_points,
                                    radius=self.fiber_rad)
        IP_1.run(inlets=ip_inlets)
        if trapping:
            IP_1.apply_trapping(outlets=outlets)

        alg_data = IP_1.plot_drainage_curve(inv_points=self.inv_points,
                                            lpf=lpf)
        plt.close('all')
        wrk.purge_object(IP_1)
        if partial:
            IP_1.return_results()
        return alg_data

    def test_apply_trapping(self):
        t = self
        t.w_inv = t.run_alg(inv_phase=t.water, def_phase=t.air)
        t.trap_data = t.run_alg(inv_phase=t.water, def_phase=t.air,
                                trapping=True)
        assert np.sum(t.trap_data[1] - t.w_inv[1]) <= 0

    def test_snap_off(self):
        t = self
        t.a_inv = t.run_alg(inv_phase=t.air, def_phase=t.water)
        t.snap_data = t.run_alg(inv_phase=t.air, def_phase=t.water,
                                snap_off=True)
        pass
        # assert np.abs(np.sum(t.snap_data[1] - t.a_inv[1])) > 0

    def test_partial(self):
        t = self
        t.imb_data = t.run_alg(inv_phase=t.air, def_phase=t.water,
                               trapping=True, partial=True)
        t.drn_data = t.run_alg(inv_phase=t.water, def_phase=t.air,
                               partial=True)
        assert t.drn_data[1][0] > 0

    def test_coop_filling(self):
        t = self
        t.w_inv = t.run_alg(inv_phase=t.water, def_phase=t.air)
        t.coop_data = t.run_alg(inv_phase=t.air, def_phase=t.water,
                                coop_fill=True)
        assert np.abs(np.sum(t.coop_data[1] - t.w_inv[1])) > 0

    def test_purcell_bi(self):
        t = self
        t.w_inv = t.run_alg(inv_phase=t.water, def_phase=t.air)
        t.bi_data = t.run_alg(inv_phase=t.water, def_phase=t.air,
                              cap_model='purcell_bi')
        assert np.abs(np.sum(t.bi_data[1] - t.w_inv[1])) > 0

    def test_sinusoidal(self):
        pass

if __name__ == '__main__':
    wrk.loglevel = 20
    t = MixedPercolationTest()
    t.setup_class()
    t.test_apply_trapping()
#    t.test_snap_off()
    t.test_partial()
    t.test_coop_filling()
    t.test_purcell_bi()
    t.test_sinusoidal()

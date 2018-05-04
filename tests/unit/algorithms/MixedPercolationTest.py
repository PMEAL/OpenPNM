import openpnm as op
import numpy as np
import scipy as sp
from openpnm.models import physics as pm
from openpnm.algorithms import MixedInvasionPercolation as mp
import matplotlib.pyplot as plt
from openpnm import topotools as tp

plt.close('all')
wrk = op.Workspace()
wrk.loglevel = 50


class MixedPercolationTest:

    def setup_class(self):
        # Create Topological Network object
        scale = 5e-5
        self.fiber_rad = 2e-6
        sp.random.seed(1)
        base_points = sp.random.random([25, 3])*scale
        self.wrk = op.core.Workspace()
        net = op.materials.VoronoiFibers(fiber_rad=2e-6,
                                         resolution=1e-6,
                                         shape=[scale, scale, scale],
                                         points=base_points,
                                         name='test')
        self.net = net
        self.prj = net.project
        self.geom = self.prj.geometries()['test_del']
        tp.trim(net, pores=net.pores('voronoi'))
        tp.trim(net, throats=net.throats()[net["throat.area"] == 0])
        tp.trim(net, throats=net.throats()[net["throat.perimeter"] == 0])
        h = net.check_network_health()
        if len(h['trim_pores']) > 0:
            tp.trim(net, pores=h['trim_pores'])
        self.air = op.phases.Air(network=net)
        self.water = op.phases.Water(network=net)
        self.water["pore.contact_angle"] = 110
        self.Swp_star = 0.25  # late pore filling
        self.air["pore.contact_angle"] = 70
        self.air["pore.surface_tension"] = self.water["pore.surface_tension"]
        self.inv_points = np.linspace(0, 30000, 31)
        # Label surfaces
        left = self.net['pore.coords'][self.net.pores('surface')][:, 0].min()
        right = self.net['pore.coords'][self.net.pores('surface')][:, 0].max()
        back = self.net['pore.coords'][self.net.pores('surface')][:, 1].min()
        front = self.net['pore.coords'][self.net.pores('surface')][:, 1].max()
        bottom = self.net['pore.coords'][self.net.pores('surface')][:, 2].min()
        top = self.net['pore.coords'][self.net.pores('surface')][:, 2].max()
        self.net['pore.left_boundary']=self.net['pore.coords'][:, 0]==left
        self.net['pore.right_boundary']=self.net['pore.coords'][:, 0]==right
        self.net['pore.back_boundary']=self.net['pore.coords'][:, 1]==back
        self.net['pore.front_boundary']=self.net['pore.coords'][:, 1]==front
        self.net['pore.bottom_boundary']=self.net['pore.coords'][:, 2]==bottom
        self.net['pore.top_boundary']=self.net['pore.coords'][:, 2]==top

    def process_physics(self, model='purcell', snap_off=True):
        prj = self.net.project
        # Clean up already made phys
        try:
            phys_air = prj.physics()['mp_phys_a']
            prj.purge_object(phys_air)
        except:
            pass
        try:
            phys_water = prj.physics()['mp_phys_w']
            prj.purge_object(phys_water)
        except:
            pass
        phys_air = op.physics.GenericPhysics(network=self.net,
                                             phase=self.air,
                                             geometry=self.geom,
                                             name='mp_phys_a')
        phys_water = op.physics.GenericPhysics(network=self.net,
                                               phase=self.water,
                                               geometry=self.geom,
                                               name='mp_phys_w')
        throat_diam = 'throat.diameter'
        pore_diam = 'pore.indiameter'
        if model == 'purcell':
            pmod = pm.capillary_pressure.purcell
            phys_water.add_model(propname='throat.capillary_pressure',
                                 model=pmod,
                                 r_toroid=self.fiber_rad,
                                 diameter=throat_diam)
            phys_air.add_model(propname='throat.capillary_pressure',
                               model=pmod,
                               r_toroid=self.fiber_rad,
                               diameter=throat_diam)
#        elif model == 'purcell_bi':
#            pmod = pm.capillary_pressure.purcell_bi
#            phys_water.add_model(propname='throat.capillary_pressure',
#                                 model=pmod,
#                                 r_toroid=self.fiber_rad,
#                                 diameter=throat_diam,
#                                 h_max=pore_diam)
#            phys_air.add_model(propname='throat.capillary_pressure',
#                               model=pmod,
#                               r_toroid=self.fiber_rad,
#                               diameter=throat_diam,
#                               h_max=pore_diam)
        elif model == 'sinusoidal':
            pmod = pm.capillary_pressure.sinusoidal
            phys_water.add_model(propname='throat.capillary_pressure',
                                 model=pmod)
            phys_air.add_model(propname='throat.capillary_pressure',
                               model=pmod)
        if snap_off:
            phys_air.add_model(propname='throat.snap_off',
                               model=pm.capillary_pressure.ransohoff_snap_off,
                               throat_diameter=throat_diam,
                               wavelength=self.fiber_rad)
            phys_air['throat.snap_off'] = np.abs(phys_air['throat.snap_off'])
        phys_air['pore.capillary_pressure'] = 0
        phys_water['pore.capillary_pressure'] = 0
        BPs = self.net.pores('surface')
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

        IP_1 = mp(network=self.net)
        IP_1.settings['partial_saturation']=partial
        IP_1.setup(phase=inv_phase,
                   def_phase=def_phase,
                   inlets=ip_inlets,
                   inlet_inv_seq=inlet_inv_seq,
                   snap_off=snap_off)
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
        self.net.project.purge_object(IP_1)
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

#    def test_purcell_bi(self):
#        t = self
#        t.w_inv = t.run_alg(inv_phase=t.water, def_phase=t.air)
#        t.bi_data = t.run_alg(inv_phase=t.water, def_phase=t.air,
#                              cap_model='purcell_bi')
#        assert np.abs(np.sum(t.bi_data[1] - t.w_inv[1])) > 0

    def test_sinusoidal(self):
        t = self
        t.purc = t.run_alg(inv_phase=t.water, def_phase=t.air)
        t.sinu = t.run_alg(inv_phase=t.water, def_phase=t.air,
                           cap_model='sinusoidal')
        assert np.abs(np.sum(t.purc[1] - t.sinu[1])) > 0

    def test_sinusoidal_coop(self):
        t = self
        t.sinu = t.run_alg(inv_phase=t.air, def_phase=t.water,
                           cap_model='sinusoidal')
        t.sinu_coop = t.run_alg(inv_phase=t.air, def_phase=t.water,
                                coop_fill=True,
                                cap_model='sinusoidal')
#        assert np.abs(np.sum(t.sinu[1] - t.sinu_coop[1])) > 0
        pass

    def test_apply_flow_rate(self):
        t = self
        pvol = np.sum(t.net['pore.volume'])
        tvol = np.sum(t.net['throat.volume'])
        tot = pvol+tvol
        t.process_physics(model='purcell', snap_off=False)
        IP_1 = mp(network=self.net)
        inlets = t.net.pores(labels=['bottom_boundary'])
        outlets = t.net.pores(labels=['top_boundary'])
        IP_1.setup(phase=t.water,
                   def_phase=t.air,
                   inlets=inlets)
        IP_1.run(outlets=outlets)
        IP_1.apply_flow(flowrate=tot)
        assert 'throat.invasion_time' in self.water.props()

if __name__ == '__main__':
    wrk.loglevel = 20
    t = MixedPercolationTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()

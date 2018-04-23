# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 11:32:32 2018

@author: Tom
"""
import openpnm as op
import numpy as np
import scipy as sp
from openpnm.models import physics as pm
from openpnm.models import geometry as gm
import matplotlib.pyplot as plt
plt.close('all')
wrk = op.core.Workspace()

cap_model = 'sinusoidal'
fiber_rad = 5e-5
snap_off = False
partial = False
coop_fill = False
trapping = False
lpf = False

sp.random.seed(0)
pn = op.network.Cubic(shape=[5, 5, 5], spacing=3e-4, name='pn11')
pn.add_boundary_pores()
proj = pn.project

Ps = pn.pores()
Ts = pn.throats()
geom = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts,
                                name='intern')
geom.add_model(propname='throat.normal', model=gm.throat_normal.pore_coords,
               regen_mode='normal')
air = op.phases.Air(network=pn)
water = op.phases.Water(network=pn)
water['throat.viscosity'] = water['pore.viscosity'][0]

mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys_air = op.physics.GenericPhysics(network=pn, phase=air,
                                     geometry=geom)
phys_water = op.physics.GenericPhysics(network=pn, phase=water,
                                       geometry=geom)
phys_water.add_model(propname='throat.conductance',
                     model=mod,
                     viscosity='throat.viscosity')

water["pore.contact_angle"] = 110
Swp_star = 0.25  # late pore filling
air["pore.contact_angle"] = 70
air["pore.surface_tension"] = water["pore.surface_tension"]
inv_points = np.linspace(0, 30000, 31)

throat_diam = 'throat.diameter'
pore_diam = 'pore.indiameter'
if cap_model == 'purcell':
    pmod = pm.capillary_pressure.purcell
elif cap_model == 'purcell_bi':
    pmod = pm.capillary_pressure.purcell_bi
elif cap_model == 'sinusoidal':
    pmod = pm.capillary_pressure.sinusoidal
phys_water.add_model(propname='throat.capillary_pressure',
                     model=pmod,
                     r_toroid=fiber_rad,
                     diameter=throat_diam,
                     h_max=pore_diam)
phys_air.add_model(propname='throat.capillary_pressure',
                   model=pmod,
                   r_toroid=fiber_rad,
                   diameter=throat_diam,
                   h_max=pore_diam)
if snap_off:
    phys_air.add_model(propname='throat.snap_off',
                       model=pm.capillary_pressure.ransohoff_snap_off,
                       diameter=throat_diam,
                       wavelength=fiber_rad)
    phys_air['throat.snap_off'] = np.abs(phys_air['throat.snap_off'])
phys_air['pore.capillary_pressure'] = 0
phys_water['pore.capillary_pressure'] = 0
inlets = pn.pores(labels=['bottom_boundary'])
outlets = pn.pores(labels=['top_boundary'])
in_step = 2
ip_inlets = [inlets[x] for x in range(0, len(inlets), in_step)]
inlet_inv_seq = -1

IP_1 = op.algorithms.MixedPercolation(network=pn)
IP_1.setup(phase=water,
           def_phase=air,
           inlets=ip_inlets,
           inlet_inv_seq=inlet_inv_seq,
           snap_off=snap_off,
           partial=partial)

if coop_fill:
    IP_1.setup_coop_filling(capillary_model=cap_model,
                            inv_points=inv_points,
                            radius=fiber_rad)
IP_1.run(inlets=ip_inlets)
if trapping:
    IP_1.apply_trapping(outlets=outlets)

alg_data = IP_1.plot_drainage_curve(inv_points=inv_points,
                                    lpf=lpf)

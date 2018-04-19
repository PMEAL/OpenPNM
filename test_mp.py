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

model = 'sinusoidal'
fiber_rad = 5e-5
snap_off = True

sp.random.seed(0)
pn = op.network.Cubic(shape=[5, 5, 5], spacing=2.5e-5, name='pn11')
pn.add_boundary_pores()
proj = pn.project

Ps = pn.pores('internal')
Ts = pn.throats('internal')
geom = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts,
                                name='intern')

Ps = pn.pores('*boundary')
Ts = pn.throats('*boundary')
boun = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts, name='bound')

pn['pore.inlets'] = pn['pore.top_boundary'].copy()
pn['pore.outlets'] = pn['pore.bottom_boundary'].copy()

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
if model == 'purcell':
    pmod = pm.capillary_pressure.purcell
elif model == 'purcell_bi':
    pmod = pm.capillary_pressure.purcell_bi
elif model == 'sinusoidal':
    pmod = pm.capillary_pressure.sinusoidal
phys_water.models.add(propname='throat.capillary_pressure',
                      model=pmod,
                      r_toroid=fiber_rad,
                      diameter=throat_diam,
                      h_max=pore_diam)
phys_air.models.add(propname='throat.capillary_pressure',
                    model=pmod,
                    r_toroid=fiber_rad,
                    diameter=throat_diam,
                    h_max=pore_diam)
if snap_off:
    phys_air.models.add(propname='throat.snap_off',
                        model=pm.capillary_pressure.ransohoff_snap_off,
                        throat_diameter=throat_diam,
                        wavelength=fiber_rad)
    phys_air['throat.snap_off'] = np.abs(phys_air['throat.snap_off'])
phys_air['pore.capillary_pressure'] = 0
phys_water['pore.capillary_pressure'] = 0
BPs = pn.pores('boundary')
NBPs = pn.find_neighbor_pores(BPs, flatten=False)
boundary_neighbors = []
for NBP in NBPs:
    boundary_neighbors.append(NBP[0])
NBPs = np.asarray(boundary_neighbors)
wPc_NBPs = phys_water["pore.capillary_pressure"][NBPs]
phys_water["pore.capillary_pressure"][BPs] = wPc_NBPs
aPc_NBPs = phys_air["pore.capillary_pressure"][NBPs]
phys_air["pore.capillary_pressure"][BPs] = aPc_NBPs
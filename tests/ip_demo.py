# -*- coding: utf-8 -*-
"""
Created on Wed May 21 16:32:22 2014

@author: MESTP
"""

import numpy as np
import OpenPNM

#==============================================================================
'''Use imported network'''
#==============================================================================

#pn = OpenPNM.Network.MatFile(name='pnMat',loglevel=30)
#pn.generate(filename='OpenPNM_net_sample_image',xtra_pore_data='type')
#geom = pn._geometries['imported']
#pn.set_pore_info(label='top',locations=pn.get_pore_indices()[pn.get_pore_data(prop='type')==6])
#pn.set_pore_info(label='bottom',locations=pn.get_pore_indices()[pn.get_pore_data(prop='type')==1])


#==============================================================================
'''Use cubic network'''
#==============================================================================

pn = OpenPNM.Network.Cubic(name='cubic_1',loglevel=30)
pn.generate(divisions=[5, 5, 5], lattice_spacing=[0.0001],add_boundaries=True)
geom = OpenPNM.Geometry.Toray090(network=pn)
geom.set_locations(pores=pn.pores('internal'),throats='all')
boun = pn.add_geometry(name='boundary_geometry',subclass='Boundary')
boun.set_locations(pores=pn.pores('boundary'))
pn.regenerate_geometries()

#==============================================================================
'''Build Fluids'''
#==============================================================================
air = OpenPNM.Fluids.Air(network=pn, name='water', loglevel=20)
air.apply_conditions(temperature=350, pressure=200000)
air.add_property(prop='electrical_conductivity',model='constant',value=5e-12)

water = OpenPNM.Fluids.Water(network=pn, name='water', loglevel=20)
water.add_property(prop='diffusivity', prop_name='DAB', model='constant', value=5e-12)

#Use Network's Fluid regeneration method
pn.regenerate_fluids()

#==============================================================================
'''Build Physics Objects'''
#==============================================================================
phys_water = OpenPNM.Physics.GenericPhysics(network=pn, fluid=water,geometry=geom,name='phys_water')
phys_water.add_property(prop='capillary_pressure', model='washburn')
#phys_water.add_property(prop='hydraulic_conductance', model='hagen_poiseuille')
#phys_water.add_property(prop='diffusive_conductance', prop_name='gdAB', model='bulk_diffusion', diffusivity='DAB')

#phys_air = OpenPNM.Physics.GenericPhysics(network=pn, fluid=air,geometry=geom, name='phys_air')
#phys_air.add_property(prop='hydraulic_conductance', model='hagen_poiseuille')
#phys_air.add_property(prop='diffusive_conductance', model='bulk_diffusion')
#phys_air.add_property(prop='electronic_conductance', model='series_resistors')

#Use Network's Physics regeneration method
pn.regenerate_physics()

#==============================================================================
'''Begin Simulations'''
#==============================================================================

ip = OpenPNM.Algorithms.InvasionPercolation(name='ip',network=pn,loglevel=30)
ip.run(inlets=pn.get_pore_indices('bottom'),outlets=pn.get_pore_indices('top'),invading_fluid=water,defending_fluid=air)
ip.update()

vis = OpenPNM.Visualization.Vtp
vis.write(pn,fluids=water)
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 08 09:43:02 2013

@author: jgostick
"""

import OpenPNM
import scipy as sp
from time import clock

import scipy.ndimage as spim
sphere = sp.ones((51,51,51),dtype=sp.bool8)
sphere[26,26,26] = 0
sphere = spim.distance_transform_edt(sphere)
template  = sphere<20

start=clock()

pn = OpenPNM.Network.GenericNetwork()

params_geo1= {'domain_size': [],  #physical network size [meters]
                'divisions': [10,10,10], #Number of pores in each direction
          'lattice_spacing': [0.1],  #spacing between pores [meters]
                'num_pores': 1000, #This is used for random networks where spacing is irrelevant
                 'template': template, #This is used for the Template based network generation
             'stats_pores' : {'name': 'weibull_min', #Each statistical package takes different params, so send as dict
                             'shape': 1.5,
                               'loc': 6e-6,
                             'scale': 2e-5},
           'stats_throats' : {'name': 'weibull_min',
                             'shape': 1.5,
                               'loc': 6e-6,
                             'scale': 2e-5},
                    'btype': [0,0,0]  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
}

#Generate Network Geometry
pn = OpenPNM.Geometry.Cubic(loglevel=40).generate(**params_geo1)
#pn = OpenPNM.Geometry.Delaunay().generate(**params)
#pn = OpenPNM.Geometry.Template().generate(**params)


#Define the fluids and set their properties
air_recipe = {       'name': 'air',
                       'Pc': 3.771e6, #Pa
                       'Tc': 132.65,  #K
                       'MW': 0.0291,  #kg/mol
              'diffusivity': {'method': 'Fuller',
                                  'MA': 31.99,
                                  'MB': 28.01,
                                  'vA': 16.6,
                                  'vB': 17.9},
                'viscosity': {'method': 'Reynolds',
                                  'uo': 0.001,
                                   'b': 0.1},
            'molar_density': {'method': 'ideal_gas',
                                   'R': 8.314},
          'surface_tension': {'method': 'na',
                               'value': [],}
}
water_recipe = {     'name': 'water',
                       'Pc': 2.206e6, #Pa
                       'Tc': 647,     #K
                       'MW': 0.0181,  #kg/mol
              'diffusivity': {'method': 'constant',
                               'value': 1e-12},
                'viscosity': {'method': 'constant',
                               'value': 0.001},
            'molar_density': {'method': 'constant',
                               'value': 44445},
          'surface_tension': {'method': 'constant',
                               'value': 0.072,}
}

#Create fluids
air = OpenPNM.Fluids.GenericFluid(loglevel=50).create(air_recipe)
water= OpenPNM.Fluids.GenericFluid(loglevel=50).create(water_recipe)

#set water and air as a fluid pair
water.set_pair(air)

#Set base conditions in the Fluids
air.pore_conditions['temperature'] = 353
air.pore_conditions['pressure'] = 101325
water.pore_conditions['temperature'] = 353
water.pore_conditions['pressure'] = 101325

water.regenerate()
air.regenerate()


OpenPNM.Physics.CapillaryPressure.set_contact_angle(water,120)
OpenPNM.Physics.CapillaryPressure.Washburn(pn,water)

OP_1 = OpenPNM.Algorithms.OrdinaryPercolation()
params_OP = { 'invading_fluid': water,
             'defending_fluid': air,
                        'npts': 50,
                   'inv_sites': [0], }
OP_1.run(pn,**params_OP)

#Assign boundary conditions
BCtypes = sp.zeros(pn.get_num_pores())
BCvalues = sp.zeros(pn.get_num_pores())
#Dirichlet
BCtypes[pn.pore_properties['type']==1] = 1
BCtypes[pn.pore_properties['type']==6] = 1
BCvalues[pn.pore_properties['type']==1] = 8e-2
BCvalues[pn.pore_properties['type']==6] = 8e-1
#Neumann
#BCtypes[pn.pore_properties['type']==1] = 1
#BCtypes[pn.pore_properties['type']==6] = 4
#BCvalues[pn.pore_properties['type']==1] = 8e-1
#BCvalues[pn.pore_properties['type']==6] = 2e-10

OpenPNM.Physics.MultiPhase.update_occupancy_OP(water,Pc=2000)
Fickian_alg = OpenPNM.Algorithms.FickianDiffusion()
Fickian_alg.set_boundary_conditions(types=BCtypes,values=BCvalues)
params_alg = {      'fluid1': air,
    'conduit_filling_method': 'strict',
                        'Pc': 0,
             }
Fickian_alg.run(pn,**params_alg)







########################################
#print "Now demonstrating the mislabeling of the Water class as 'air'"
#pn = OpenPNM.Geometry.Cubic(loglevel=40).generate(**params_geom1)
#air = OpenPNM.Fluids.Water().create(name='air')
#water= OpenPNM.Fluids.GenericFluid(loglevel=50).create(water_recipe)
#air.assign_to_network(pn)
#water.assign_to_network(pn)
#print ''
#print 'current pore conditions:'
#for i in pn.pore_conditions.keys():
#    print i,'=',pn.pore_conditions[i]
#print ''
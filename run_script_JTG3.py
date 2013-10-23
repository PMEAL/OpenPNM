# -*- coding: utf-8 -*-
"""
Created on Fri Mar 08 09:43:02 2013

@author: jgostick
"""

import OpenPNM
import scipy as sp
from time import clock
import copy

import scipy.ndimage as spim
sphere = sp.ones((51,51,51),dtype=sp.bool8)
sphere[26,26,26] = 0
sphere = spim.distance_transform_edt(sphere)
template  = sphere<20

start=clock()

pn = OpenPNM.Network.GenericNetwork()

params_geom1 =     {'domain_size': [],  #physical network size [meters]
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
pn = OpenPNM.Geometry.Cubic(loglevel=40).generate(**params_geom1)
#pn = OpenPNM.Geometry.Delaunay().generate(**params)
#pn = OpenPNM.Geometry.Template().generate(**params)

#Set Base Conditions in the Network
pn.pore_conditions['temperature'] = 353
pn.pore_conditions['pressure'] = 101325

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
}
#Create fluids
air = OpenPNM.Fluids.GenericFluid(loglevel=50).create(air_recipe)
water= OpenPNM.Fluids.GenericFluid(loglevel=50).create(water_recipe)

#Assign fluids to network
air.assign_to_network(pn)
water.assign_to_network(pn)
print ''
print 'current pore conditions:'
for i in pn.pore_conditions.keys():
    print i,'=',pn.pore_conditions[i]
print ''
#Run some algorithms that change base conditions
#blah, blah, blah...
print 'changing temp and pressure...'
pn.pore_conditions['temperature'] = 333
pn.pore_conditions['pressure'] = 201325

#Update fluids
air.refresh_in_network(pn)
water.refresh_in_network(pn)
print ''
print 'current pore conditions:'
for i in pn.pore_conditions.keys():
    print i,'=',pn.pore_conditions[i]
print ''

print "Swapping out fluid 'water' with a similar 'solution'"
#Create a new fluid that is similar to water
solution_recipe = pn.phases['water']
#Subtly change something about it
solution_recipe['viscosity']['value'] = 0.0015
solution = OpenPNM.Fluids.GenericFluid().create(solution_recipe)
solution.rename_fluid('solution')
#Add this fluid to the network
solution.assign_to_network(pn)
#Remove water
water.remove_from_network(pn)

print ''
print 'current pore conditions:'
for i in pn.pore_conditions.keys():
    print i,'=',pn.pore_conditions[i]
print ''

##################################################################
#Now, try to do everything above without creating a fluid object
print 'Working from the network exclusively'
#Generate Network Geometry
pn = OpenPNM.Geometry.Cubic(loglevel=40).generate(**params_geom1)
pn.pore_conditions['temperature'] = 353
pn.pore_conditions['pressure'] = 101325

#Assign fluids to network
pn.assign_fluid(air_recipe)
pn.assign_fluid(water_recipe)

print ''
print 'current pore conditions:'
for i in pn.pore_conditions.keys():
    print i,'=',pn.pore_conditions[i]
print ''
#Run some algorithms that change base conditions
#blah, blah, blah...
print 'changing temp and pressure...'
pn.pore_conditions['temperature'] = 333
pn.pore_conditions['pressure'] = 201325

#Update fluids
pn.refresh_all_fluids()

print ''
print 'current pore conditions:'
for i in pn.pore_conditions.keys():
    print i,'=',pn.pore_conditions[i]
print ''

print "Swapping out fluid 'water' with a similar 'solution'"
#Create a new fluid that is similar to water
solution_recipe = copy.deepcopy(pn.phases['water'])
#Subtly change something about it
solution_recipe['viscosity']['value'] = 0.0015
solution_recipe['name'] = 'solution'
#Add this fluid to the network
pn.assign_fluid(solution_recipe)
#Remove water
pn.remove_fluid('water')

print ''
print 'current pore conditions:'
for i in pn.pore_conditions.keys():
    print i,'=',pn.pore_conditions[i]
print ''
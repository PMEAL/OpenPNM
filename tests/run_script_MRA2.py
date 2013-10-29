# -*- coding: utf-8 -*-
"""
Created on Tue Jul 09 10:54:37 2013

@author: Mahmoudreza Aghighi
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

#Set Base Conditions in the Network
pn.pore_conditions['temperature'] = 353
pn.pore_conditions['pressure'] = 101325

#Define the fluids and set their properties
params_air = {       'name': 'air',
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
             'wettability' : 'wp',
            'molar_density': {'method': 'ideal_gas',
                                   'R': 8.413},
}
params_water = {     'name': 'water',
                       'Pc': 2.206e6, #Pa
                       'Tc': 647,     #K
                       'MW': 0.0181,  #kg/mol
              'diffusivity': {'method': 'constant',
                               'value': 1e-12},
                'viscosity': {'method': 'constant',
                               'value': 0.001},
             'wettability' : 'wp',
            'molar_density': {'method': 'constant',
                               'value': 44445},
}
#Create fluids
air = OpenPNM.Fluids.GenericFluid(loglevel=50).create(params_air)
water= OpenPNM.Fluids.GenericFluid().create(params_water)

#Assign boundary conditions
BCtypes = sp.zeros(pn.get_num_pores())
BCvalues = sp.zeros(pn.get_num_pores())
#Dirichlet
BCtypes[pn.pore_properties['type']==1] = 1
BCtypes[pn.pore_properties['type']==6] = 1
BCvalues[pn.pore_properties['type']==1] = 8e-2
BCvalues[pn.pore_properties['type']==6] = 8e-1
##Neumann
#BCtypes[pn.pore_properties['type']==1] = 1
#BCtypes[pn.pore_properties['type']==6] = 4
#BCvalues[pn.pore_properties['type']==1] = 8e-1
#BCvalues[pn.pore_properties['type']==6] = -2e-10

Fickian_alg = OpenPNM.Algorithms.FickianDiffusion()
Fickian_alg.set_boundary_conditions(types=BCtypes,values=BCvalues)
params_alg = {
                     'fluid':air,
    'conduit_filling_method': 'strict',
                        'Pc': 0,
             }
Fickian_alg.run(pn,**params_alg)



##Permeability
#BCtypes[pn.pore_properties['type']==1] = 1
#BCtypes[pn.pore_properties['type']==6] = 1
#BCvalues[pn.pore_properties['type']==1] = 1e6
#BCvalues[pn.pore_properties['type']==6] = 3e6
#
#permeability_alg = OpenPNM.Algorithms.Permeability()
#permeability_alg.set_boundary_conditions(types=BCtypes,values=BCvalues)
#params_alg = {
#                'fluid_name': 'water',
#    'conduit_filling_method': 'strict',
#                        'Pc': 0,
#             }
#permeability_alg.run(pn,**params_alg)

print clock()-start,"seconds."
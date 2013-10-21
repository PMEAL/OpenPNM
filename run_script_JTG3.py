# -*- coding: utf-8 -*-
"""
Created on Fri Mar 08 09:43:02 2013

@author: jgostick
"""

import OpenPNM
import scipy as sp
from time import clock
import pprint

import scipy.ndimage as spim
sphere = sp.ones((51,51,51),dtype=sp.bool8)
sphere[26,26,26] = 0
sphere = spim.distance_transform_edt(sphere)
template  = sphere<20

start=clock()

pn = OpenPNM.Network.GenericNetwork()

params =     {'domain_size': [],  #physical network size [meters]
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
pn = OpenPNM.Geometry.Cubic(loglevel=50).generate(**params)
#pn = OpenPNM.Geometry.Delaunay().generate(**params)
#pn = OpenPNM.Geometry.Template().generate(**params)

#Set Base Conditions in the Network
pn.pore_conditions['temperature'] = 353
pn.pore_conditions['pressure'] = 101325

#Define the fluids and set their properties
params_air = {       'name': 'air',
                    'phase': 1,
                     'type': 'compressible_gas',
                      'Pc' : 3.771e6, #Pa
                      'Tc' : 132.65,  #K
                      'MW' : 0.0291, #kg/mol
              'diffusivity': {'method': 'Fuller',
                                  'MA': 31.99,
                                  'MB': 28.01,
                                  'vA': 16.6,
                                  'vB': 17.9},
                'viscosity': {'method': 'Reynolds',
                                  'uo': 0.001,
                                   'b': 0.1},
            'molar_density': {'method': 'ideal_gas',
                                   'R': 8.413},
}
params_water = {     'name': 'water',
                    'phase': 0,
                     'type': 'incompressible_liquid',
                      'Pc' : 2.2064e6, #Pa
                      'Tc' : 647,      #K
                      'MW' : 0.0181,   #kg/mol
              'diffusivity': {'method': 'constant',
                                 'DAB': 1e-12},
                'viscosity': {'method': 'constant',
                                  'mu': 0.001},
            'molar_density': {'method': 'constant',
                                   'c': 44445},
}

air = OpenPNM.Fluids.GenericFluid().create(**params_air)
water = OpenPNM.Fluids.GenericFluid().create(**params_water)

#Now associate fluids with a network







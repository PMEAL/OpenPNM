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

params =     {'domain_size': [],  #physical network size [meters]
                'divisions': [10,10,10], #Number of pores in each direction
          'lattice_spacing': [0.1],  #spacing between pores [meters]
                'num_pores': 1000, #This is used for random networks where spacing is irrelevant
                 'template': template, #This is used for the Template based network generation
'stats_pores' :     {'name': 'weibull_min', #Each statistical package takes different params, so send as dict
                    'shape': 1.5,
                      'loc': 6e-6,
                    'scale': 2e-5},
'stats_throats' :   {'name': 'weibull_min',
                    'shape': 1.5,
                      'loc': 6e-6,
                    'scale': 2e-5},
                    'btype': [0,0,0]  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
}

#Generate Network Geometry
pn = OpenPNM.Geometry.Cubic(loglevel=50).generate(**params)
#pn = OpenPNM.Geometry.Delaunay().generate(**params)
#pn = OpenPNM.Geometry.Template().generate(**params)

#Define the Fluids and set their properties
air_dict = {    'name':  'air',
               'phase':  1,
         'diffusivity':  2.09e-5,
           'viscosity':  1.73e-5,
        'molardensity':  40.9}
water_dict = {  'name':  'water',
               'phase':  2,
         'diffusivity':  1.0e-20,
           'viscosity':  1.0e-3,
        'molardensity':  5.56e4,
         }
air = OpenPNM.Fluids.GenericFluid(air_dict)
air.update(pn)
water = OpenPNM.Fluids.GenericFluid(water_dict)
water.update(pn)
water.add_property('bob',5)

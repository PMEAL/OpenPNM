# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 15:35:56 2013

@author: debarungupta
"""

import OpenPNM
import scipy as sp
from time import clock
import scipy.ndimage as spim


params = {
'domain_size'           : [0.001,0.001,0.0004],  #physical network size [meters]
'divisions'             : [], #Number of pores in each direction
'lattice_spacing'       : [.0001],  #spacing between pores [meters]
'num_pores'             : 1000, #This is used for random networks where spacing is irrelevant
'template'              : template, #This is used for the Template based network generation
'psd_info'   : {'name'  : 'weibull_min', #Each statistical package takes different params, so send as dict
                'shape' : 1.5,
                'loc'   : 6e-6,
                'scale' : 2e-5},
'tsd_info'   : {'name'  : 'weibull_min',
                'shape' : 1.5,
                'loc'   : 6e-6,
                'scale' : 2e-5},
'btype'                 : [0,0,0],  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
}

pn = OpenPNM.Geometry.Cubic().generate(**params)

'''
boundaries = sp.array(sp.repeat(params,6)) # Creates a list of 
temp = params['divisions']
number_of_sides = 6
mult = -1*(sp.eye(number_of_sides/2)-1) # This creates the compliment of an eye matrix

for i in range(3):
    boundary_dims = (temp*mult[i % 3,:]) + sp.eye(3)[i % 3,:]
    boundaries[i]['divisions'] = boundary_dims.tolist()
    print boundaries[i]['divisions']

OpenPNM.Generators.Cubic(loglevel=10,**params).generate()
'''
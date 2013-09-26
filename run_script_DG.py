# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 15:35:56 2013

@author: debarungupta
"""

import OpenPNM
import scipy as sp
from time import clock

params = {
#'domain_size': [0.001,0.001,0.0004],  #physical network size [meters]
'divisions': [20,15,10], #Number of pores in each direction
'lattice_spacing': 1.0,  #spacing between pores [meters]
#'num_pores': 1000, #This is used for random networks where spacing is irrelevant
'btype': [0,0,0],  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
#The following parameters are for the statistical functions used to generate
#pore and throat size distributions.  See scipy.stats for options and parameters.
'psd_dist': 'weibull_min',  #distribution to be used
'psd_shape': 1.5,  #spread or skew of distribution
'psd_loc': 6e-6,  #minimum pore diameter [meters]
'psd_scale': 2e-5,  #controls maximum pore size
'tsd_dist': 'weibull_min',  #distribution to be used
'tsd_shape': 1.5,  #spread or skew of distribution
'tsd_loc': 6e-6,  #minimum throat diameter [meters]
'tsd_scale': 2e-5  #controls maximum throat size
}

start=clock()
#pn = OpenPNM.Generators.Cubic(loglevel=10,**params).generate()

boundaries = sp.array(sp.repeat(params,6)) # Creates a list of 
temp = params['divisions']
number_of_sides = 6
mult = -1*(sp.eye(number_of_sides/2)-1) # This creates the compliment of an eye matrix

for i in range(3):
    boundary_dims = (temp*mult[i % 3,:]) + sp.eye(3)[i % 3,:]
    boundaries[i]['divisions'] = boundary_dims.tolist()
    print boundaries[i]['divisions']

OpenPNM.Generators.Cubic(loglevel=10,**params).generate()
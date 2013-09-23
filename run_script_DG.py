# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 21:13:50 2013

@author: debarungupta
"""


import OpenPNM
import numpy as np
from time import clock

params = {
#'domain_size': [0.001,0.001,0.0004],  #physical network size [meters]
'divisions': [20,20,20], #Number of pores in each direction
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
# Test Comment
start=clock()
pn1 = OpenPNM.GEN.Cubic(loglevel=10,**params).generate()
pn2 = OpenPNM.GEN.Cubic(loglevel=10,**params).generate()

# throat_properties: Volume
# throat_properties: Diameter
# throat_properties: Numbering
# throat_properties: Generate connections
size_of_sparse_1 = np.max(pn1.throat_properties['connections'])
full_sparse = np.zeros(np.repeat(size_of_sparse_1,2))   # We are attempting to convert the ijv into a binary sparse to make the connections.

# throat_properties: Length
# throat_properties: Seed
# throat_properties: Type
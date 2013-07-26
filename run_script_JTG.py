# -*- coding: utf-8 -*-
"""
Created on Fri Mar 08 09:43:02 2013

@author: jgostick
"""

import OpenPNM
import numpy as np
import matplotlib.pyplot as plt
from time import clock

params = {
'domain_size': [0.001,0.001,0.0004],  #physical network size [meters]
#'divisions': [10,10,10], #Number of pores in each direction
'lattice_spacing': 0.00005,  #spacing between pores [meters]
'num_pores': 1000, #This is used for random networks where spacing is irrelevant
'btype': [1,1,0],  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
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
pn = OpenPNM.GEN.Cubic(loglevel=20,**params).generate()
#pn = OpenPNM.GEN.Delaunay(loglevel=10,**params).generate()

pn.throat_properties['Pc_entry'] = -4*0.072*np.cos(np.radians(105))/pn.throat_properties['diameter']  #This should be set somewhere else
exp1 = OpenPNM.ALG.OrdinaryPercolationAlgorithm(pn, npts=100, inv_faces=[1]).run()
pn.update()

exp2= OpenPNM.ALG.FicksLawDiffusion(pn.exp1,a)

#Write network to vtk file for visualization in Paraview
OpenPNM.IO.NetToVtp(pn,'C:/Users/jeff/Dropbox/Flash Sync/Code/OpenPNM/30PythonPackage/trunk/OpenPNM/IO/test.vtk')

print clock()-start,"seconds."

vis = OpenPNM.VIS.Vis2D()
vis.overview(pn)
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 08 09:43:02 2013

@author: jgostick
"""

import OpenPNM
import numpy as np
import matplotlib.pyplot as plt
from time import clock

import scipy.ndimage as spim
img = spim.imread('C:\Users\jeff\Dropbox\Flash Sync\Code\Git\OpenPNM\snippets\CL.jpg')
img = img[0:100,0:100,0]<30

params = {
#'domain_size': [0.001,0.001,0.0004],  #physical network size [meters]
'divisions': [60,60,60], #Number of pores in each direction
'lattice_spacing': 1,  #spacing between pores [meters]
#'num_pores': 1000, #This is used for random networks where spacing is irrelevant
'btype': [1,1,0],  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
'template': img,
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
pn = OpenPNM.Geometry.Template(loglevel=10,**params).generate()

#pn.throat_properties['Pc_entry'] = -4*0.072*np.cos(np.radians(105))/pn.throat_properties['diameter']  #This should be set somewhere else
#inlets = [0]
#outlets = pn.get_num_pores()
#exp1 = OpenPNM.ALG.InvasionPercolationAlgorithm(pn, loglevel = 10, npts=100, inlets=inlets, outlets=outlets).run()
#exp2 = OpenPNM.ALG.OrdinaryPercolationAlgorithm(pn, npts=50, inv_sites=inlets).run()
#pn.update()

#Write network to vtk file for visualization in Paraview
#import os
#OpenPNM.IO.NetToVtp(pn,os.path.abspath(os.path.dirname(__file__))+'\OpenPNM\\IO\\test.vtk')

print clock()-start,"seconds."

#vis = OpenPNM.VIS.Vis2D()
#vis.overview(pn)

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 05 08:55:04 2013
"""

import OpenPNM
import numpy as np
import matplotlib.pyplot as plt
from time import clock

import scipy.ndimage as spim
img = np.random.rand(60,60,60)>0.0005
img = spim.distance_transform_edt(img)<=4

import scipy.ndimage as spim
tmp = spim.imread('C:\Users\Jeff\Pictures\Picture2.tif')
tmp = tmp[35:105,:,3]==0
#plt.imshow(tmp)
img = np.kron(np.ones((3,1,1)),tmp) #Make image deeper
img = np.transpose(img,axes=[1,2,0])
#img = tmp

params = {
#'domain_size': [0.001,0.001,0.0004],  #physical network size [meters]
'divisions': [60,60,60], #Number of pores in each direction
'lattice_spacing': 1,  #spacing between pores [meters]
#'num_pores': 1000, #This is used for random networks where spacing is irrelevant
'image_shape': img,
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
gn = OpenPNM.GEN.Custom(loglevel=10,**params)
gn.generate()
pn = gn.get_net()
#gn.generate_pore_property_from_image(np.array((imgsld),dtype=np.int8),'rand')

neighborPs = pn.get_neighbor_pores(np.r_[0:pn.get_num_pores()],flatten=False)
tmp = np.zeros((pn.get_num_pores(),))
for i in np.r_[0:pn.get_num_pores()]:
    tmp[i] = np.size(neighborPs[i])
pn.pore_properties['core_shell']=np.array(tmp<6,dtype=np.int8)*1 + np.array(tmp==6,dtype=np.int8)*2

#pn = OpenPNM.GEN.Cubic(loglevel=20,**params).generate()
#pn = OpenPNM.GEN.Delaunay(loglevel=10,**params).generate()

pn.throat_properties['Pc_entry'] = OpenPNM.PHYS.CapillaryPressure.Washburn(pn,sigma,theta)
pn.throat_properties['Pc_entry']= -4*0.072*np.cos(np.radians(105))/pn.throat_properties['diameter']  #This should be set somewhere else
inlets = [200]
#outlets = [pn.get_num_pores()-1]
#exp1 = OpenPNM.ALG.InvasionPercolationAlgorithm(pn, loglevel = 10, npts=100, inlets=inlets, outlets=outlets).run()
#exp2 = OpenPNM.ALG.OrdinaryPercolationAlgorithm(pn, npts=10, inv_sites=inlets).run()
#pn.update()

#Write network to vtk file for visualization in Paraview
#import os
#OpenPNM.IO.NetToVtp(pn,os.path.abspath(os.path.dirname(__file__))+'\OpenPNM\\IO\\test.vtk')

print clock()-start,"seconds."

#vis = OpenPNM.VIS.Vis2D()
#vis.overview(pn)
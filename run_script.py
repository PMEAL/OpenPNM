# -*- coding: utf-8 -*-
"""
Created on Fri Mar 08 09:43:02 2013

@author: jgostick
"""

import OpenPNM
from time import clock
import scipy as sp

import scipy.ndimage as spim
sphere = sp.ones((51,51,51),dtype=sp.bool8)
sphere[26,26,26] = 0
sphere = spim.distance_transform_edt(sphere)
template = sphere<20

params = {
'domain_size'               : [1,1,1],  #physical network size [meters]
'divisions'                 : [10,10,10], #Number of pores in each direction
'lattice_spacing'           : [],  #spacing between pores [meters]
'num_pores'                 : 1000, #This is used for random networks where spacing is irrelevant
'template'                  : template, #This is used for the Template based network generation
'stats_pores'   : {  'name' : 'weibull_min', #Each statistical package takes different params, so send as dict
                    'shape' : 1.5,
                      'loc' : 6e-6,
                    'scale' : 2e-5},
'stats_throats' : {  'name' : 'weibull_min',
                    'shape' : 1.5,
                      'loc' : 6e-6,
                    'scale' : 2e-5},
'btype'                     : [0,0,0],  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
}

start=clock()
pn = OpenPNM.Geometry.Cubic(loglevel=10).generate(**params)
#pn = OpenPNM.Geometry.Delaunay().generate(**params)
#pn = OpenPNM.Geometry.Template().generate(**params)

pn.throat_properties['Pc_entry'] = OpenPNM.Physics.CapillaryPressure.Washburn(pn,0.072,110)
inlets = [0]
#exp1 = OpenPNM.Algorithms.InvasionPercolation(pn, loglevel=10, npts=100, inlets=inlets, outlets=outlets).run()
exp2 = OpenPNM.Algorithms.OrdinaryPercolation().run(pn, npts=50, inv_sites=inlets)

#Write network to vtk file for visualization in Paraview
#import os
#filename = os.path.abspath(os.path.dirname(__file__))+'\LocalFiles\\test.vtk'
#OpenPNM.Visualization.VTK().write(pn,filename)

#print pn
print clock(),"seconds."

#vis = OpenPNM.Algorithms.Vis2D()
#vis.overview(pn)

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 08 09:43:02 2013

@author: jgostick
"""

import OpenPNM
import scipy as sp
import matplotlib.pyplot as plt
from time import clock

params = {
'domain_size'       : [0.001,0.001,0.0004],  #physical network size [meters]
'divisions'         : [], #Number of pores in each direction
'lattice_spacing'   : [.0001],  #spacing between pores [meters]
'num_pores'         : 1000, #This is used for random networks where spacing is irrelevant
'psd_info'          : {'name'  : 'weibull_min',
                       'shape' : 1.5,
                       'loc'   : 6e-6,
                       'scale' : 2e-5},
'tsd_info'          : {'name'  : 'weibull_min',
                       'shape' : 1.5,
                       'loc'   : 6e-6,
                       'scale' : 2e-5},
'btype'             : [0,0,0],  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
}

start=clock()
pn = OpenPNM.Geometry.Cubic().generate(**params)

#pn = OpenPNM.Geometry.Delaunay(loglevel=10,**params).generate()

pn.throat_properties['Pc_entry'] = OpenPNM.Physics.CapillaryPressure().Washburn(pn,0.072,110)
inlets = [0]
#exp1 = OpenPNM.Algorithms.InvasionPercolation(pn, loglevel=10, npts=100, inlets=inlets, outlets=outlets).run()
exp2 = OpenPNM.Algorithms.OrdinaryPercolation(pn, npts=50, inv_sites=inlets).run()
pn.update()

#Write network to vtk file for visualization in Paraview
#import os
#OpenPNM.Visualization.NetToVtp(pn,os.path.abspath(os.path.dirname(__file__))+'\OpenPNM\\IO\\test.vtk')

print clock()-start,"seconds."

#vis = OpenPNM.Algorithms.Vis2D()
#vis.overview(pn)

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 15:35:56 2013

@author: debarungupta
"""

import OpenPNM
import scipy as sp
from time import clock
import scipy.ndimage as spim
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


# Parameters unique to all matricies.
Nx = 3
Ny = 3
Nz = 3

network_main = {
'psd_info'   : {'name'  : 'weibull_min', #Each statistical package takes different params, so send as dict
                'shape' : 1.5,
                'loc'   : 6e-6,
                'scale' : 2e-5},
'tsd_info'   : {'name'  : 'weibull_min',
                'shape' : 1.5,
                'loc'   : 6e-6,
                'scale' : 2e-5},
'btype'                 : [0,0,0],  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
'lattice_spacing'       : [.01],  #spacing between pores [meters]
'divisions'             : [Nx,Ny,Nz]
}

# Parameters specific to individual matricies.

#Generate the main pore networks.
pn1 = OpenPNM.Geometry.Cubic().generate(**network_main)
pn2 = OpenPNM.Geometry.Cubic().generate(**network_main)

#Add boundaries to the networks
OpenPNM.Geometry.Cubic()._generate_boundaries(pn1,**network_main)

#Stitch the networks
OpenPNM.Geometry.Cubic().stitch_network(pn1,pn2,stitch_side = 'bottom') # can be stitched to top, bottom, left right etc.

OpenPNM.Geometry.GenericGeometry().plot_net(pn1)
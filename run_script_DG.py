# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 15:35:56 2013

@author: debarungupta
"""

import OpenPNM
import scipy as sp
from time import clock
import scipy.ndimage as spim

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

#Generate the main pore network.
pn = OpenPNM.Geometry.Cubic().generate(**network_main)
OpenPNM.Geometry.Cubic()._generate_boundaries(pn,**network_main)
# Call add boundaries(**params). Add boundares keeps calls generate and creates new networks based on the Nx, Ny, Nz.
    #Add boundaries then calls stitch after generating. Stitch is put into cubic.stitch(). 

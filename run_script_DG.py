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
params = {
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
}

Nx = 4
Ny = 4
Nz = 4
# Parameters specific to individual matricies.
main_params = {'divisions' : [Nx,Ny,Nz]}

#Generate the main pore network.
network_main = dict(params.items() + main_params.items())
pn = OpenPNM.Geometry.Cubic().generate(**network_main)

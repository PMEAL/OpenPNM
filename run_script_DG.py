# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 15:35:56 2013

@author: debarungupta
"""

import OpenPNM
import scipy as sp
from time import clock
import scipy.ndimage as spim


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
}

specific_params_1 = {
'divisions'             : [10,10,10], #Number of pores in each direction
'lattice_spacing'       : [.01],  #spacing between pores [meters]
}

specific_params_2 = {
'divisions'             : [20,20,20], #Number of pores in each direction
'lattice_spacing'       : [.01],  #spacing between pores [meters]
}

network1_props = dict(params.items() + specific_params_1.items())
network2_props = dict(params.items() + specific_params_2.items())

pn1 = OpenPNM.Geometry.Cubic().generate(**network1_props)
pn2 = OpenPNM.Geometry.Cubic().generate(**network2_props)

OpenPNM.Geometry.GenericGeometry().translate_coordinates(pn2,displacement = [0,0,pn1.pore_properties['coords'][:,2].max()])
OpenPNM.Geometry.GenericGeometry().stitch(pn1,pn2)


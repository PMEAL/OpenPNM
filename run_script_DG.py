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
'lattice_spacing'       : [.01],  #spacing between pores [meters]
}

main_params = {'divisions' : [4,4,4]}
side_1 = {'divisions' : [2,4,4]}
side_2 = {'divisions' : [4,1,4]}
side_3 = {'divisions' : [4,4,1]}
side_4 = {'divisions' : [1,4,4]}
side_5 = {'divisions' : [4,1,4]}
side_6 = {'divisions' : [4,4,1]}

network_main = dict(params.items() + main_params.items())
pn = OpenPNM.Geometry.Cubic().generate(**network_main)

net_side1 = dict(params.items() + side_1.items())
net_side2 = dict(params.items() + side_2.items())
net_side3 = dict(params.items() + side_3.items())
net_side4 = dict(params.items() + side_4.items())
net_side5 = dict(params.items() + side_5.items())
net_side6 = dict(params.items() + side_6.items())

pn1 = OpenPNM.Geometry.Cubic().generate(**net_side1)
pn2 = OpenPNM.Geometry.Cubic().generate(**net_side2)
pn3 = OpenPNM.Geometry.Cubic().generate(**net_side3)
pn4 = OpenPNM.Geometry.Cubic().generate(**net_side4)
pn5 = OpenPNM.Geometry.Cubic().generate(**net_side5)
pn6 = OpenPNM.Geometry.Cubic().generate(**net_side6)

OpenPNM.Geometry.GenericGeometry().translate_coordinates(pn1,displacement = [pn.pore_properties['coords'][:,0].max(),0,0])
OpenPNM.Geometry.GenericGeometry().translate_coordinates(pn2,displacement = [0,pn.pore_properties['coords'][:,1].max(),0])
OpenPNM.Geometry.GenericGeometry().translate_coordinates(pn3,displacement = [0,0,pn.pore_properties['coords'][:,2].max()])
OpenPNM.Geometry.GenericGeometry().translate_coordinates(pn4,displacement = [-1*pn.pore_properties['coords'][:,0].min(),0,0])
OpenPNM.Geometry.GenericGeometry().translate_coordinates(pn5,displacement = [0,-1*pn.pore_properties['coords'][:,1].min(),0])
OpenPNM.Geometry.GenericGeometry().translate_coordinates(pn6,displacement = [0,0,-1*pn.pore_properties['coords'][:,2].min()])

OpenPNM.Geometry.GenericGeometry().stitch(pn,pn1,edge=1)
OpenPNM.Geometry.GenericGeometry().stitch(pn,pn2,edge=2)
OpenPNM.Geometry.GenericGeometry().stitch(pn,pn3,edge=3)
OpenPNM.Geometry.GenericGeometry().stitch(pn,pn4,edge=4)
OpenPNM.Geometry.GenericGeometry().stitch(pn,pn5,edge=5)
OpenPNM.Geometry.GenericGeometry().stitch(pn,pn6,edge=6)

# Can we add extra keys in the dictionary to return the input parmeters so we can access it for stitching?


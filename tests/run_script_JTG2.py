# -*- coding: utf-8 -*-
"""
Created on Fri Mar 08 09:43:02 2013

@author: jgostick
"""

import OpenPNM
import scipy as sp
from time import clock
import pprint

import scipy.ndimage as spim
sphere = sp.ones((51,51,51),dtype=sp.bool8)
sphere[26,26,26] = 0
sphere = spim.distance_transform_edt(sphere)
#template  = sphere<20

start=clock()

params =     {'domain_size': [],  #physical network size [meters]
                'divisions': [10,10,10], #Number of pores in each direction
          'lattice_spacing': [0.1],  #spacing between pores [meters]
                'num_pores': 1000, #This is used for random networks where spacing is irrelevant
                 'template': template, #This is used for the Template based network generation
'stats_pores' :     {'name': 'weibull_min', #Each statistical package takes different params, so send as dict
                    'shape': 1.5,
                      'loc': 6e-6,
                    'scale': 2e-5},
'stats_throats' :   {'name': 'weibull_min',
                    'shape': 1.5,
                      'loc': 6e-6,
                    'scale': 2e-5},
                    'btype': [0,0,0]  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
}

#Generate Network Geometry
pn = OpenPNM.Geometry.Cubic(loglevel=50).generate(**params)
#pn = OpenPNM.Geometry.Delaunay().generate(**params)
#pn = OpenPNM.Geometry.Template().generate(**params)

#Define the Fluids and set their properties
air = {'name': 'air'}
OpenPNM.Fluids.Diffusivity.set_as(air,2.09e-5)
OpenPNM.Fluids.Viscosity.set_as(air,1.73e-5)
OpenPNM.Fluids.MolarDensity.set_as(air,40.90)
water = {'name': 'water'}
water['critical_temperature'] = 500
water['critical_pressure'] = 1e10
OpenPNM.Fluids.Diffusivity.set_as(water,1.0e-20)
OpenPNM.Fluids.Viscosity.set_as(water,1.0e-3)
OpenPNM.Fluids.MolarDensity.set_as(water,5.56e4)
OpenPNM.Fluids.SurfaceTension.set_as(water,air,0.072)
solid = {'name': 'solid'}
OpenPNM.Fluids.SurfaceTension.set_as(solid,air,0.01)
OpenPNM.Fluids.SurfaceTension.set_as(solid,water,0.02)

#Apply Pore Scale Physics
OpenPNM.Physics.MassTransport.DiffusiveConductance(pn,air)
OpenPNM.Physics.MassTransport.DiffusiveConductance(pn,water)
OpenPNM.Physics.CapillaryPressure.set_contact_angle(water,120)
OpenPNM.Physics.CapillaryPressure.Washburn(pn,water,air)

pprint.pprint(air,   width=30, depth=1)
pprint.pprint(water, width=30, depth=1)

#Perform Algorithms
inlets = sp.r_[0:pn.get_num_pores()]
mask = pn.pore_properties['type']==2
inlets = pn.pore_properties['numbering'][mask]
OpenPNM.Algorithms.OrdinaryPercolation(loglevel=10).run(net=pn, invading_fluid=water, defending_fluid=air, npts=50, inv_sites=inlets)

#OpenPNM.Algorithms.InvasionPercolation(loglevel=10).run(pn,inlets=[0],outlets=[1],end_condition='breakthrough',timing='ON',report=20)

#Write network to vtk file for visualization in Paraview
#import os
#fname = os.path.abspath(os.path.dirname(__file__))+'\LocalFiles\\test.vtk'
#OpenPNM.Visualization.VTK(loglevel=50).write(pn,fname)
#OpenPNM.Visualization.Plots.Capillary_Pressure_Curve(pn)

print clock()-start,"seconds."


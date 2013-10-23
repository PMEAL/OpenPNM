# -*- coding: utf-8 -*-
"""
Created on Tue Jul 09 10:54:37 2013

@author: Mahmoudreza Aghighi
"""
import OpenPNM
import scipy as sp
from time import clock
import copy

import scipy.ndimage as spim
sphere = sp.ones((51,51,51),dtype=sp.bool8)
sphere[26,26,26] = 0
sphere = spim.distance_transform_edt(sphere)
template  = sphere<20

start=clock()

pn = OpenPNM.Network.GenericNetwork()

params_geom1 =     {'domain_size': [],  #physical network size [meters]
                'divisions': [10,10,10], #Number of pores in each direction
          'lattice_spacing': [0.1],  #spacing between pores [meters]
                'num_pores': 1000, #This is used for random networks where spacing is irrelevant
                 'template': template, #This is used for the Template based network generation
             'stats_pores' : {'name': 'weibull_min', #Each statistical package takes different params, so send as dict
                             'shape': 1.5,
                               'loc': 6e-6,
                             'scale': 2e-5},
           'stats_throats' : {'name': 'weibull_min',
                             'shape': 1.5,
                               'loc': 6e-6,
                             'scale': 2e-5},
                    'btype': [0,0,0]  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
}

#Generate Network Geometry
pn = OpenPNM.Geometry.Cubic(loglevel=40).generate(**params_geom1)
#pn = OpenPNM.Geometry.Delaunay().generate(**params)
#pn = OpenPNM.Geometry.Template().generate(**params)

#Set Base Conditions in the Network
pn.pore_conditions['temperature'] = 353
pn.pore_conditions['pressure'] = 101325

#Define the fluids and set their properties
params_air = {       'name': 'air',
                       'Pc': 3.771e6, #Pa
                       'Tc': 132.65,  #K
                       'MW': 0.0291,  #kg/mol
              'diffusivity': {'method': 'Fuller',
                                  'MA': 31.99,
                                  'MB': 28.01,
                                  'vA': 16.6,
                                  'vB': 17.9},
                'viscosity': {'method': 'Reynolds',
                                  'uo': 0.001,
                                   'b': 0.1},
             'wettability' : 'wp',
            'molar_density': {'method': 'ideal_gas',
                                   'R': 8.413},
}
params_water = {     'name': 'water',
                       'Pc': 2.206e6, #Pa
                       'Tc': 647,     #K
                       'MW': 0.0181,  #kg/mol
              'diffusivity': {'method': 'constant',
                               'value': 1e-12},
                'viscosity': {'method': 'constant',
                               'value': 0.001},
             'wettability' : 'wp',
            'molar_density': {'method': 'constant',
                               'value': 44445},
}
#Create fluids
air = OpenPNM.Fluids.GenericFluid(loglevel=50).create(params_air)
water= OpenPNM.Fluids.GenericFluid().create(params_water)

#Assign fluids to network
air.assign_to_network(pn)
water.assign_to_network(pn)

#inlets = sp.nonzero(pn.pore_properties['type']==1)[0]
#pn.throat_properties['Pc_entry'] = OpenPNM.Physics.CapillaryPressure().Washburn(pn,0.072,110)  #This should be set somewhere else
#OP = OpenPNM.Algorithms.OrdinaryPercolation(pn, npts=50, inv_sites=inlets)
#OP.run()

BCtypes = sp.zeros(pn.get_num_pores())
BCvalues = sp.zeros(pn.get_num_pores())
BCtypes[pn.pore_properties['type']==1]=1
BCtypes[pn.pore_properties['type']==6]=1
BCvalues[pn.pore_properties['type']==1]=0.5
BCvalues[pn.pore_properties['type']==6]=0.2

Alg1=OpenPNM.Algorithms.FickianDiffusion()
Alg1.set_boundary_conditions(types=BCtypes,values=BCvalues)
Alg1.run(pn,fluid_name='air')

#pn.update()
#setattr(pn,"Total_Conc",C)
#setattr(pn,"Diff_Coefficient",D)
#setattr(pn,"divisions",params['divisions'])
##setattr(pn,"lattice_spacing",lattice_spacing)
#
#
#"FOR IP"
##IP = OpenPNM.Algorithms.InvasionPercolationAlgorithmTiming(net=pn,loglevel=20,end_condition='breakthrough')
##IP.run()
##pn.set_pore_property('INVADED PORES',np.multiply(np.ones(pn.get_num_pores()),pn.pore_properties['IP_Pseq']<max(pn.pore_properties['IP_Pseq']/2)))
##MT = OpenPNM.Algorithms.FickianDiffusion(pn,Alg='IP',Psequence=[max(pn.pore_properties['IP_Pseq'])/2])
#
#"FOR OP"
#pn.pore_properties['Pc_invaded'] = sp.random.rand(pn.get_num_pores())
#P = 0
##pn.pore_properties['INVADED PORES'] = sp.multiply(sp.ones(pn.get_num_pores()),pn.pore_properties['Pc_invaded']<P)
#MT = OpenPNM.Algorithms.FickianDiffusion(pn)
#
#"FOR None"
##list1=range(401,492)
##list2=range(501,592)
##list3=range(601,692)
##LIST=list1+list2+list3
##Pores =np.array(LIST)
##pn.set_pore_property('INVADED PORES',np.multiply(np.ones(pn.get_num_pores()),np.in1d(np.array(range(pn.get_num_pores())),Pores)))
##MT = OpenPNM.Algorithms.FickianDiffusion(pn,Alg='None',Pinvaded=Pores)
#BCtypes = sp.zeros(self._net.get_num_pores())
#BCvalues = sp.zeros(self._net.get_num_pores())
#
#
#MT.run()
##Write network to vtk file for visualization in Paraview

print clock()-start,"seconds."
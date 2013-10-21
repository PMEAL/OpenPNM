# -*- coding: utf-8 -*-
"""
Created on Tue Jul 09 10:54:37 2013

@author: Mahmoudreza Aghighi
"""

import OpenPNM
import scipy as sp
from time import clock

r"""
Except Network generator and Algorithms, all of the following parameters and variables 
should be in the other dictionaries or classes.
In here, just for the sake of simplicity, they have been set in the run-script file
"""

C = 0.9
D = 2.09e-9
params = {
'domain_size'               : [1,1,1],  #physical network size [meters]
'divisions'                 : [10,10,10], #Number of pores in each direction
'lattice_spacing'           : [],  #spacing between pores [meters]
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
pn = OpenPNM.Geometry.Cubic().generate(**params)
inlets = sp.nonzero(pn.pore_properties['type']==1)[0]
pn.throat_properties['Pc_entry'] = OpenPNM.Physics.CapillaryPressure().Washburn(pn,0.072,110)  #This should be set somewhere else
OP = OpenPNM.Algorithms.OrdinaryPercolation(pn, npts=50, inv_sites=inlets)
OP.run()
pn.update()
setattr(pn,"Total_Conc",C)
setattr(pn,"Diff_Coefficient",D)
setattr(pn,"divisions",params['divisions'])
#setattr(pn,"lattice_spacing",lattice_spacing)


"FOR IP"
#IP = OpenPNM.Algorithms.InvasionPercolationAlgorithmTiming(net=pn,loglevel=20,end_condition='breakthrough')
#IP.run()
#pn.set_pore_property('INVADED PORES',np.multiply(np.ones(pn.get_num_pores()),pn.pore_properties['IP_Pseq']<max(pn.pore_properties['IP_Pseq']/2)))
#MT = OpenPNM.Algorithms.FickianDiffusion(pn,Alg='IP',Psequence=[max(pn.pore_properties['IP_Pseq'])/2])

"FOR OP"
pn.pore_properties['Pc_invaded'] = sp.random.rand(pn.get_num_pores())
P = 0
#pn.pore_properties['INVADED PORES'] = sp.multiply(sp.ones(pn.get_num_pores()),pn.pore_properties['Pc_invaded']<P)
MT = OpenPNM.Algorithms.FickianDiffusion(pn)

"FOR None"
#list1=range(401,492)
#list2=range(501,592)
#list3=range(601,692)
#LIST=list1+list2+list3
#Pores =np.array(LIST)
#pn.set_pore_property('INVADED PORES',np.multiply(np.ones(pn.get_num_pores()),np.in1d(np.array(range(pn.get_num_pores())),Pores)))
#MT = OpenPNM.Algorithms.FickianDiffusion(pn,Alg='None',Pinvaded=Pores)



MT.run()
#Write network to vtk file for visualization in Paraview

print clock()-start,"seconds."
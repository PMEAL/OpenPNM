# -*- coding: utf-8 -*-
"""
Created on Tue Jul 09 10:54:37 2013

@author: Mahmoudreza Aghighi
"""

import OpenPNM
import numpy as np
from time import clock

start=clock()
r"""
Except Network generator and Algorithms, all of the following parameters and variables 
should be in the other dictionaries or classes.
In here, just for the sake of simplicity, they have been set in the run-script file
"""

C = 0.9
D = 2.09e-9
divisions=[10,10,10]
lattice_spacing=.0005
pn = OpenPNM.Geometry.Cubic(divisions=[10,10,10],lattice_spacing=.0005,loglevel=20,btype = [0,0,0]).generate()
setattr(pn,"Total_Conc",C)
setattr(pn,"Diff_Coefficient",D)
setattr(pn,"divisions",divisions)
setattr(pn,"lattice_spacing",lattice_spacing)
pn.throat_properties['Pc_entry'] = -4*0.072*np.cos(np.radians(105))/pn.throat_properties['diameter']  #This should be set somewhere else

"FOR IP"
#IP = OpenPNM.Algorithms.InvasionPercolationAlgorithmTiming(net=pn,loglevel=20,end_condition='breakthrough')
#IP.run()
#pn.set_pore_property('INVADED PORES',np.multiply(np.ones(pn.get_num_pores()),pn.pore_properties['IP_Pseq']<max(pn.pore_properties['IP_Pseq']/2)))
#MT = OpenPNM.Algorithms.FickianDiffusion(pn,Alg='IP',Psequence=[max(pn.pore_properties['IP_Pseq'])/2])

"FOR OP"
OP = OpenPNM.Algorithms.OrdinaryPercolationAlgorithm(pn, npts=100, inv_faces=[1])
OP.run()
P = 0.04
pn.set_pore_property('INVADED PORES',np.multiply(np.ones(pn.get_num_pores()),pn.pore_properties['Pc_invaded']<P))
MT = OpenPNM.Algorithms.FickianDiffusion(pn,Alg='OP',Pressure=[P])

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
OpenPNM.Visualization.NetToVtp(net = pn,filename="MassTransfer.vtp")

print clock()-start,"seconds."
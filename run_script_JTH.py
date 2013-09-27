# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 00:00:01 2013
╛
@author: jhinebaugh
"""

import OpenPNM
from time import clock
import scipy as sp

from time import clock


start=clock()
print "="*50
print "= Example: Create random network and run an invasion\n= percolation algorithm"
print "-"*50
print "- * generate a simple cubic network"    
#sp.random.seed(1)
pn = OpenPNM.Generators.Cubic(divisions=[10,10,10],lattice_spacing=.0005,loglevel=20,btype = [0,0,0]).generate()
#pn = OpenPNM.Generators.Delaunay(domain_size=[30,30,10],num_pores = 5000 ,btype = [1,1,0]).generate()
print "+"*50
print "Sample generated at t =",clock()-start,"seconds."
print "+"*50
#print "- * Assign pore volumes"
#pore_volumes=sp.random.rand(pn.get_num_pores())
#pore_volumes[range(pn.get_num_pores([0]),pn.get_num_pores())]=0
#pn.set_pore_property(name='volume',ndarray=pore_volumes,columns = None)  
print '- * Assign boundary pore volumes = 0'
pn.pore_properties['diameter'][pn.pore_properties['type']>0] = 0
    

print "- * Define inlet and outlet faces"
inlets = sp.nonzero(pn.pore_properties['type']==1)[0]
outlets = sp.nonzero(pn.pore_properties['type']==6)[0]
inlets2 = sp.unique((inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
                   inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
                   inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
                   inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
                   inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))]))
print inlets2
#print "- * assign random pore and throat diameters"
#pn.pore_properties['diameter'] = sp.random.rand(pn.get_num_pores(),1)
#pn.throat_properties['diameter'] = sp.random.rand(pn.get_num_throats(),1)

print "- * Run Invasion percolation algorithm"
IP = OpenPNM.Algorithms.InvasionPercolationAlgorithmTiming(net=pn,inlets=inlets2[1],outlets=outlets,loglevel=20,loggername="TestInvPercAlg")
IP.run()
print "+"*50
print "IP completed at t =",clock()-start,"seconds."
print "+"*50
print "- * Save output to IP.vtp"
OpenPNM.Visualization.NetToVtp(net = pn,filename="IP.vtp")

print "="*50
print "Program Finished at t = ",clock()-start,"seconds."
print "="*50

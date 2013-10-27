  # -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 00:00:01 2013
╛
@author: jhinebaugh
"""

import OpenPNM
from time import clock
import scipy as sp
clock()
print 'Creating 1,000 pore cubic network'
cubic = OpenPNM.Geometry.Cubic().generate(domain_size=[10,10,10],lattice_spacing=[1.0])
print cubic
print 'Finished at time =',sp.round_(clock(),2),'seconds'
print ''
print 'Creating 100 pore Delauny network'
delaunay = OpenPNM.Geometry.Delaunay().generate(domain_size=[10,10,10],num_pores=100)
print delaunay
print 'Finished at time =',sp.round_(clock(),2),'seconds'
print ''
print 'Creating #### pore imported network'
matfile = OpenPNM.Geometry.MatFile().generate(filename='large_network')
print matfile
print 'Finished at time =',sp.round_(clock(),2),'seconds'
print ''

print "Adding 'air' and 'water' for each network type"
air_cubic = OpenPNM.Fluids.Air().create(fluid_name='air_cubic')
water_cubic = OpenPNM.Fluids.Water().create(fluid_name='water_cubic')
air_delaunay = OpenPNM.Fluids.Air().create(fluid_name='air_delaunay')
water_delaunay = OpenPNM.Fluids.Water().create(fluid_name='water_delaunay')
air_matfile = OpenPNM.Fluids.Air().create(fluid_name='air_matfile')
water_matfile = OpenPNM.Fluids.Water().create(fluid_name='water_matfile')

print ''
print "pairing fluids"
air_cubic.set_pair(water_cubic)
air_delaunay.set_pair(water_delaunay)
air_matfile.set_pair(water_matfile)
print ''
print 'Running IP algorithm for cubic'
inlets = cubic.pore_properties['numbering'][cubic.pore_properties['type']==1]
outlets = cubic.pore_properties['numbering'][cubic.pore_properties['type']==6]
OpenPNM.Algorithms.InvasionPercolation().run(cubic,invading_fluid=water_cubic,inlets=inlets,outlets=outlets)
print ''
print 'Running IP algorithm for delaunay'
inlets = delaunay.pore_properties['numbering'][delaunay.pore_properties['type']==1]
outlets = delaunay.pore_properties['numbering'][delaunay.pore_properties['type']==6]
OpenPNM.Algorithms.InvasionPercolation().run(delaunay,invading_fluid=water_delaunay,inlets=inlets,outlets=outlets)
print ''
print 'Running IP algorithm for matfile'
inlets = matfile.pore_properties['numbering'][matfile.pore_properties['type']==1]
outlets = matfile.pore_properties['numbering'][matfile.pore_properties['type']==6]
OpenPNM.Algorithms.InvasionPercolation().run(matfile,invading_fluid=water_matfile,inlets=inlets,outlets=outlets)
print ''
print 'Running OP algorithm for cubic'
inlets = cubic.pore_properties['numbering'][cubic.pore_properties['type']==1]
outlets = cubic.pore_properties['numbering'][cubic.pore_properties['type']==6]
OpenPNM.Algorithms.OrdinaryPercolation().run(water_cubic,air_cubic,inlets=inlets,outlets=outlets,num_points=25)
print ''
print 'Running OP algorithm for delaunay'
inlets = delaunay.pore_properties['numbering'][delaunay.pore_properties['type']==1]
outlets = delaunay.pore_properties['numbering'][delaunay.pore_properties['type']==6]
OpenPNM.Algorithms.OrdinaryPercolation().run(water_delaunay,air_delaunay,inlets=inlets,outlets=outlets,num_points=25)
print ''
print 'Running OP algorithm for matfile'
inlets = matfile.pore_properties['numbering'][matfile.pore_properties['type']==1]
outlets = matfile.pore_properties['numbering'][matfile.pore_properties['type']==6]
OpenPNM.Algorithms.OrdinaryPercolation().run(water_matfile,air_matfile,inlets=inlets,outlets=outlets,num_points=25)




#
#
#print "="*50
#print "= Example: Create random network and run an invasion\n= percolation algorithm"
#print "-"*50
#print "- * generate a simple cubic network"    
##sp.random.seed(1)
#pn = OpenPNM.Geometry.Cubic(divisions=[10,10,10],lattice_spacing=.0005,loglevel=20,btype = [0,0,0]).generate()
##pn = OpenPNM.Geometry.Delaunay(domain_size=[30,30,10],num_pores = 5000 ,btype = [1,1,0]).generate()
#print "+"*50
#print "Sample generated at t =",clock()-start,"seconds."
#print "+"*50
##print "- * Assign pore volumes"
##pore_volumes=sp.random.rand(pn.get_num_pores())
##pore_volumes[range(pn.get_num_pores([0]),pn.get_num_pores())]=0
##pn.set_pore_property(name='volume',ndarray=pore_volumes,columns = None)  
#print '- * Assign boundary pore volumes = 0'
#pn.pore_properties['diameter'][pn.pore_properties['type']>0] = 0
#    
#
#print "- * Define inlet and outlet faces"
#inlets = sp.nonzero(pn.pore_properties['type']==1)[0]
#outlets = sp.nonzero(pn.pore_properties['type']==6)[0]
#inlets2 = sp.unique((inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
#                   inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
#                   inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
#                   inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))],
#                   inlets[sp.random.randint(sp.size(inlets,0))],inlets[sp.random.randint(sp.size(inlets,0))]))
#print inlets2
##print "- * assign random pore and throat diameters"
##pn.pore_properties['diameter'] = sp.random.rand(pn.get_num_pores(),1)
##pn.throat_properties['diameter'] = sp.random.rand(pn.get_num_throats(),1)
#
#print "- * Run Invasion percolation algorithm"
#IP = OpenPNM.Algorithms.InvasionPercolation(net=pn,inlets=inlets2[1],outlets=outlets,loglevel=20,loggername="TestInvPercAlg")
#IP.run()
#print "+"*50
#print "IP completed at t =",clock()-start,"seconds."
#print "+"*50
#print "- * Save output to IP.vtp"
#OpenPNM.Visualization.NetToVtp(net = pn,filename="IP.vtp")
#
#print "="*50
#print "Program Finished at t = ",clock()-start,"seconds."
#print "="*50

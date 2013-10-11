#! /usr/bin/env python

# Author: OpenPNM Development Team
# License: TBD
# Copyright (c) 2012, TBD
r"""

"""

import OpenPNM
import scipy as sp
from time import clock
print "="*50
print "= Example: Create random network and run an invasion\n= percolation algorithm"
print "-"*50
print "- * generate a simple cubic network"
pn = OpenPNM.Generators.Cubic().generate()

start=clock()

print "- * Define inlet and outlet faces"
inlets = sp.nonzero(pn.pore_properties['type']==1)[0]
outlets = sp.nonzero(pn.pore_properties['type']==6)[0]

print "- * assign random pore and throat diameters"
pn.pore_properties['dia'] = sp.rand(pn.get_num_pores(),1)
pn.throat_properties['dia'] = sp.rand(pn.get_num_throats(),1)

print "- * Run Invasion percolation algorithm"
IP = OpenPNM.Algorithms.InvasionPercolation(net=pn,inlets=[inlets],outlets=outlets,loggername="InvPercAlg",loglevel=30)
IP.run()

print "- * Save output to IP.vtp"
OpenPNM.Visualization.NetToVtp(net = pn,filename="IP.vtp")

print "="*50
print "Delta T: ",clock()-start,"seconds."
print "="*50
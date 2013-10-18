# -*- coding: utf-8 -*-

# Getting started

import OpenPNM
gn1 = OpenPNM.Geometry.Cubic()
pn1 = gn1.generate()
print pn1
pn2 = OpenPNM.Geometry.Cubic().generate()

#Network Architecture
pn = OpenPNM.Geometry.Cubic().generate(divisions=[3,3,3],lattice_spacing=[1])
print pn
pn.print_overview()
print pn.pore_properties.keys()
pn.pore_conditions['temperature'] = 80.0
pn.pore_conditions['temperature'] = 80.1
pn.pore_conditions['pressure'] = 101325
gas_constant = 8.314
pn.pore_conditions['molar_density'] = pn.pore_conditions['pressure']/gas_constant/pn.pore_conditions['temperature']
pn.print_overview()
import scipy as sp
z_mean = sp.mean(pn.pore_properties['coords'][:,2])
mask = pn.pore_properties['coords'][:,2] < z_mean
low_pores = pn.pore_properties['numbering'][mask]
print low_pores
print pn.pore_properties['type'][low_pores]
Nt = pn.get_num_throats()
values = sp.random.rand(Nt,)*4 + 1 # 1 < ratios < 5
pn.throat_properties['aspect_ratio'] = values
print 'There are',pn.get_num_pores(),'pores in the network.'
print 'There are',pn.get_num_throats(),'throats in the network.'
print 'Pore 5 has the following neighbors:',pn.get_neighbor_pores(5)
print 'Pore 5 has the following throats:',pn.get_neighbor_throats(5)
print 'Pore 5 has',pn.get_num_neighbors(5),'neighbors.'
print 'Throat 6 connects pore',pn.get_connected_pores(6)[0],'to pore',pn.get_connected_pores(6)[1],'.'
print 'Throat',pn.get_connecting_throat(0,1),'connects pore 0 to pore 1.'


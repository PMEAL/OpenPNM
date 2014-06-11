# -*- coding: utf-8 -*-
"""
Drainange Curve on a cubic network example

@author: pmtgt
"""
import OpenPNM
import matplotlib.pyplot as plt
" Create Network Object "
" Simple Cubic "
#pn = OpenPNM.Network.Cubic(name='test').generate(lattice_spacing=[0.0001],divisions=[10,10,10],add_boundaries=True)
" Cyclindrical "
#pn = OpenPNM.Network.Cylinder(name='test').generate(radius=0.001, length=0.001, lattice_spacing=0.0001,add_boundaries=True)
" Spherical "
#pn = OpenPNM.Network.Sphere(name='test').generate(radius=0.001, lattice_spacing=0.0001,add_boundaries=True)
" Delaunay Network "
pn = OpenPNM.Network.Delaunay(name='del')
pn.generate(num_pores=250, domain_size=[0.0001,0.0001,0.0001],add_boundaries=True) # Good pore stats 
#pn.generate(num_pores=1000, domain_size=[0.0002,0.0002,0.0001],add_boundaries=True) # Good pore stats 
#pn.generate(num_pores=4000, domain_size=[0.0004,0.0004,0.0001],add_boundaries=True) # Same as above on larger domain
#loc = pn.pores()
" Create Geometry Object and invoke "
#geo = OpenPNM.Geometry.Stick_and_Ball(network=pn,name='basic')
" Initialis and build the geometry object "
#geom = OpenPNM.Geometry.GenericGeometry(network=pn)  # instantiate geometry object
geom = OpenPNM.Geometry.Voronoi(network=pn,loglevel=20,name='vor')  # instantiate geometry object
" Add desired methods to geometry object "# All put inside new class
geom.set_locations(pores=pn.pores('internal'),throats='all')
pn.regenerate_geometries()
" Create Fluid Objects and invoke "
air = OpenPNM.Fluids.Air(network=pn,loglevel=20,name='air')
air.regenerate()
water = OpenPNM.Fluids.Water(network=pn,loglevel=20,name='water')
water.add_property(prop='diffusivity',prop_name='DAB',model='constant',value=5e-12)
water.regenerate()
#==============================================================================
'''Build Physics Objects'''
#==============================================================================
phys_water = OpenPNM.Physics.GenericPhysics(network=pn, fluid=water,geometry=geom)
phys_water.add_property(prop='capillary_pressure', model='washburn')
phys_water.add_property(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_water.add_property(prop='diffusive_conductance', model='bulk_diffusion', diffusivity='DAB')

phys_air = OpenPNM.Physics.GenericPhysics(network=pn, fluid=air,geometry=geom)
phys_air.add_property(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_air.add_property(prop='diffusive_conductance', model='bulk_diffusion')
#phys_air.add_property(prop='electronic_conductance', model='series_resistors')

#Use Network's Physics regeneration method
pn.regenerate_physics()
" Change contact angle "
#water.add_method(prop='contact_angle',model='constant',value=100)
#water.regenerate()
" Run a drainage simulation "
" Create and algorithm object, define injection sites and run "

OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(loglevel=20,network=pn)
a = pn.pores(labels=['bottom','boundary'],mode='intersection')
OP_1.setup(invading_fluid=water,defending_fluid=air,inlets=a,npts=100)
OP_1.run()

#OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn,name='OP_1')
#injection_sites = pn.get_pore_indices(labels='bottom')
#OP_1.setup(invading_fluid='water',defending_fluid='air',inlets=injection_sites,npts=20)
#OP_1.run()
" Plot Results "
OP_1.plot_drainage_curve()
" Get Results for a certain Capillary Pressure "
OP_1.update(Pc=20000)
" Output "
OpenPNM.Visualization.Vtp.write(filename='test.vtp',fluids=[air,water,pn],network=pn)

pore_vols = pn.get_pore_data(prop='volume')
dom_vol = pn._Lx*pn._Ly*pn._Lz
porosity = sum(pore_vols)/dom_vol
print("Porosity: "+str(porosity))
throat_d = pn.get_throat_data(prop='diameter')
pore_d = pn.get_pore_data(prop='diameter')
" Plot histograms of pore stats"
plt.figure()
plt.title("Histogram of Pore Volumes")
num,bins,patches = plt.hist(pore_vols,bins=20,normed=1,histtype='bar',color='r')
plt.show()

plt.figure()
plt.title("Histogram of Pore Diameters")
num,bins,patches = plt.hist(pore_d,bins=20,normed=1,histtype='bar',color='b')
plt.show()

plt.figure()
plt.title("Histogram of Throat Diameters")
num,bins,patches = plt.hist(throat_d,bins=20,normed=1,histtype='bar',color='g')
plt.show()
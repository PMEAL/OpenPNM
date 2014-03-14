import OpenPNM
import scipy as sp

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.MatFile(name='pnMat',loglevel=10).generate(filename='standard_cubic_5x5x5.mat')
#pn = OpenPNM.Network.Cubic(name='cubic_1',loglevel=10).generate(divisions=[10, 10, 10], lattice_spacing=[0.0001],add_boundaries=True)
#pn = OpenPNM.Network.Delaunay(name='random_1',loglevel=10).generate(num_pores=1500,domain_size=[100,100,30])
#pn = OpenPNM.Network.Template(name='template_1',loglevel=10).generate(template=sp.ones((4,4),dtype=int),lattice_spacing=0.001)
#pn = OpenPNM.Network.Sphere(name='sphere_1',loglevel=10).generate(radius=5,lattice_spacing=1)
#pn = OpenPNM.Network.Cylinder(name='cylinder1',loglevel=10).generate(radius=10,length=5,lattice_spacing=1)
#pn = OpenPNM.Network.TestNet()

#==============================================================================
'''Build Geometry'''
#==============================================================================
GDL_geom = pn._geometry[0]
GDL_geom.add_method(prop='throat_length',model='constant',value=1)
GDL_geom.regenerate()
#GDL_pores = sp.r_[0:1500]
#GDL_throats = pn.find_neighbor_throats(GDL_pores,mode='intersection')
#GDL_geom = OpenPNM.Geometry.Stick_and_Ball(network=pn, name='GDL', pnums=GDL_pores, tnums=GDL_throats)
#GDL_geom.regenerate()
#------------------------------------------------------------------------------
#MPL_pores = sp.r_[1500:pn.num_pores()]
#MPL_throats = pn.find_neighbor_throats(MPL_pores,mode='intersection')
#MPL_geom = OpenPNM.Geometry.Stick_and_Ball(network=pn, name='MPL', pnums=MPL_pores, tnums=MPL_throats)
#MPL_geom.regenerate()
#------------------------------------------------------------------------------
#interface_throats = pn.find_interface_throats(['GDL','MPL'])
#interface_geom = OpenPNM.Geometry.GenericGeometry(network=pn, name='interface',tnums=interface_throats)
#interface_geom.add_method(prop='throat_seed',model='neighbor_min')
#interface_geom.add_method(prop='throat_diameter',model='cylinder',name='weibull_min',shape=2.5,loc=6e-6,scale=2e-5)
#interface_geom.add_method(prop='throat_length',model='straight')
#interface_geom.add_method(prop='throat_volume',model='cylinder')
#interface_geom.add_method(prop='throat_vector',model='pore_to_pore')
#interface_geom.add_method(prop='throat_area',model='cylinder')
#interface_geom.add_method(prop='throat_surface_area',model='cylinder')
#interface_geom.regenerate()

#==============================================================================
'''Build Fluids'''
#==============================================================================
air = OpenPNM.Fluids.Air(network=pn, loglevel=10,init_cond={'temperature':300, 'pressure':100000})
air.apply_ICs(init_cond={'temperature':350, 'pressure':200000})  # experimental feature
air.regenerate()

water = OpenPNM.Fluids.Water(network=pn,loglevel=10)
water.add_method(prop='diffusivity',prop_name='DAB',model='constant',value=5e-12)
water.regenerate()
#
##==============================================================================
#'''Build Physics Objects'''
##==============================================================================
phys_water_GDL = OpenPNM.Physics.GenericPhysics(network=pn, fluid='water',geometry=GDL_geom,name='phys_water_GDL')
phys_water_GDL.add_method(prop='capillary_pressure', model='purcell', r_toroid=1e-5)
phys_water_GDL.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_water_GDL.add_method(prop='diffusive_conductance', prop_name='gdAB', model='bulk_diffusion', diffusivity='DAB')
phys_water_GDL.regenerate()

phys_air_GDL = OpenPNM.Physics.GenericPhysics(network=pn, fluid=air,geometry=GDL_geom, name='phys_air_GDL')
phys_air_GDL.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_air_GDL.add_method(prop='diffusive_conductance', model='bulk_diffusion')
phys_air_GDL.regenerate()
## ----------------------------------------------------------------------------------------------
#phys_water_MPL = OpenPNM.Physics.GenericPhysics(network=pn, fluid=water,geometry='MPL',name='phys_water_MPL')
#phys_water_MPL.add_method(prop='capillary_pressure', model='purcell', r_toroid=1e-5)
#phys_water_MPL.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
#phys_water_MPL.add_method(prop='diffusive_conductance', prop_name='gdAB', model='bulk_diffusion', diffusivity='DAB')
#phys_water_MPL.regenerate()
#
#phys_air_MPL = OpenPNM.Physics.GenericPhysics(network=pn, fluid=air, geometry='MPL',name='phys_air_MPL')
#phys_air_MPL.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
#phys_air_MPL.add_method(prop='diffusive_conductance', model='bulk_diffusion')
#phys_air_MPL.regenerate()
## ----------------------------------------------------------------------------------------------
#phys_water_interface = OpenPNM.Physics.GenericPhysics(network=pn, fluid=water,geometry='interface',name='phys_water_interface')
#phys_water_interface.add_method(prop='capillary_pressure', model='purcell', r_toroid=1e-5)
#phys_water_interface.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
#phys_water_interface.add_method(prop='diffusive_conductance', prop_name='gdAB', model='bulk_diffusion', diffusivity='DAB')
#phys_water_interface.regenerate()
#
#phys_air_interface = OpenPNM.Physics.GenericPhysics(network=pn, fluid=air, geometry='interface',name='phys_air_interface')
#phys_air_interface.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
#phys_air_interface.add_method(prop='diffusive_conductance', model='bulk_diffusion')
#phys_air_interface.regenerate()

#==============================================================================
'''Begin Simulations'''
#==============================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#------------------------------------------------------------------------------
#Initialize algorithm object
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(loglevel=10,loggername='OP',name='OP_1',network=pn)
a = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
OP_1.setup(invading_fluid='water',defending_fluid='air',inlets=a,npts=20)
OP_1.run()


#b = pn.get_pore_indices(labels=['top','boundary'],mode='intersection')
#OP_1.evaluate_trapping(outlets=b)
#OP_1.plot_drainage_curve()

##-----------------------------------------------------------------------------
#'''Perform an Injection Experiment (InvasionPercolation)'''
##-----------------------------------------------------------------------------
##Initialize algorithm object
#IP_1 = OpenPNM.Algorithms.InvasionPercolation(loglevel=10,name='IP_1',network=pn)
#face = pn.get_pore_indices(labels=['right','boundary'],mode='intersection',indices=False)
#quarter = sp.rand(pn.num_pores(),)<.1
#inlets = pn.get_pore_indices()[face&quarter]
#outlets = pn.get_pore_indices('left')
#IP_1.run(invading_fluid=water,defending_fluid=air,inlets=inlets,outlets=outlets)

##----------------------------------------------------------------------
#'''Perform Fickian Diffusion'''
##----------------------------------------------------------------------
## Updating data based on the result of Percolation Algorithms
OP_1.update(Pc=3000)
#IP_1.update()
###----------------------------------------------------------------------
### Initializing diffusion algorithm
Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loglevel=10, loggername='Fickian', name='Fickian_alg',network=pn)
###----------------------------------------------------------------------------
### Assign Dirichlet boundary conditions to some of the surface pores
##BC1
BC1_pores = pn.get_pore_indices(labels=['top','boundary'],mode='intersection')
Fickian_alg.set_pore_info(label='Dirichlet', locations=BC1_pores)
BC1_values = 0.6
Fickian_alg.set_pore_data(prop='BCval', data=BC1_values, locations=BC1_pores)
## BC2
BC2_pores = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
Fickian_alg.set_pore_info(label='Dirichlet', locations=BC2_pores)
BC2_values = 0.2
Fickian_alg.set_pore_data(prop='BCval', data=BC2_values, locations=BC2_pores)
###----------------------------------------------------------------------------
### Assign Neumann and Dirichlet boundary conditions to some of the surface pores
### BC1
#BC1_pores = pn.get_pore_indices(labels=['top','boundary'],mode='intersection')
#Fickian_alg.set_pore_info(label='Dirichlet',locations=BC1_pores)
#BC1_values = 0.5
#Fickian_alg.set_pore_data(prop='BCval',data=BC1_values,locations=BC1_pores)
### BC2
#BC2_pores = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
#Fickian_alg.set_pore_info(label='Neumann_rate_group',locations=BC2_pores)
#BC2_values = 2e-9
#Fickian_alg.set_pore_data(prop='BCval',data=BC2_values,locations=BC2_pores)
###----------------------------------------------------------------------------
### Assign Dirichlet boundary condition to some of the surface pores and 
### Neumann boundary condition to all of the internal pores(individually) 
###BC0
#BC0_pores = pn.get_pore_indices()[-sp.in1d(pn.get_pore_indices(),pn.get_pore_indices(['top','bottom']))]
#Fickian_alg.set_pore_info(label='Neumann_rate_single',locations=BC0_pores)
#BC0_values = 7e-12
#Fickian_alg.set_pore_data(prop='BCval',data=BC0_values,locations=BC0_pores)
###BC1
#BC1_pores = pn.get_pore_indices(labels=['top','boundary'],mode='intersection')
#Fickian_alg.set_pore_info(label='Dirichlet',locations=BC1_pores)
#BC1_values = 0.6
#Fickian_alg.set_pore_data(prop='BCval',data=BC1_values,locations=BC1_pores)
### BC2
#BC2_pores = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
#Fickian_alg.set_pore_info(label='Dirichlet',locations=BC2_pores)
#BC2_values = 0.2
#Fickian_alg.set_pore_data(prop='BCval',data=BC2_values,locations=BC2_pores)
###----------------------------------------------------------------------------
### Assign Dirichlet boundary condition to some of the surface pores and 
### Neumann boundary condition to some of the internal pores(to the cluster not individually)
###BC0
#BC0_pores = [500,501,502,503,504]
#Fickian_alg.set_pore_info(label='Neumann_rate_group',locations=BC0_pores)
#BC0_values = 5e-7
#Fickian_alg.set_pore_data(prop='BCval',data=BC0_values,locations=BC0_pores)
###BC1
#BC1_pores = pn.get_pore_indices(labels=['top','boundary'],mode='intersection')
#Fickian_alg.set_pore_info(label='Dirichlet',locations=BC1_pores)
#BC1_values = 0.4
#Fickian_alg.set_pore_data(prop='BCval',data=BC1_values,locations=BC1_pores)
### BC2
#BC2_pores = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
#Fickian_alg.set_pore_info(label='Dirichlet',locations=BC2_pores)
#BC2_values = 0.3
#Fickian_alg.set_pore_data(prop='BCval',data=BC2_values,locations=BC2_pores)
###----------------------------------------------------------------------------
### Assign Dirichlet boundary condition to some of the surface pores and 
### Neumann insulated boundary condition to some of the internal pores
###BC0
#BC0_pores = sp.r_[500:530]
#Fickian_alg.set_pore_info(label='Neumann_insulated',locations=BC0_pores)
###BC1
#BC1_pores = pn.get_pore_indices(labels=['top','boundary'],mode='intersection')
#Fickian_alg.set_pore_info(label='Dirichlet',locations=BC1_pores)
#BC1_values = 0.4
#Fickian_alg.set_pore_data(prop='BCval',data=BC1_values,locations=BC1_pores)
### BC2
#BC2_pores = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
#Fickian_alg.set_pore_info(label='Dirichlet',locations=BC2_pores)
#BC2_values = 0.1
#Fickian_alg.set_pore_data(prop='BCval',data=BC2_values,locations=BC2_pores)
###----------------------------------------------------------------------------
### Run simulation
Fickian_alg.run(active_fluid=air)
Fickian_alg.update()
###-----------------------------------------------------------------------
###Export to VTK
OpenPNM.Visualization.VTK().write(net=pn, fluids=[air,water])
### Capillary pressure curve and Overview histograms
#OpenPNM.Visualization.__Plots__.Capillary_Pressure_Curve(net=pn,fluid=water)
#OpenPNM.Visualization.__Plots__.Overview(net=pn)
import OpenPNM

#==============================================================================
'''Build Topological Network'''
#==============================================================================
#pn = OpenPNM.Network.Cubic(name='cubic_1').generate(divisions=[35,35,35],lattice_spacing=[0.0001])
pn = OpenPNM.Network.Delaunay(name='random_1',loglevel=10).generate(num_pores=100,domain_size=[100,100,100])
#pn = OpenPNM.Network.TestNet()

#==============================================================================
'''Build Geometry'''
#==============================================================================
geom = OpenPNM.Geometry.Stick_and_Ball(network=pn,name='stick_and_ball',locations=pn.get_pore_indices())
geom.regenerate()

#==============================================================================
'''Build Fluids'''
#==============================================================================
air = OpenPNM.Fluids.Air(network=pn)
air.regenerate()

water = OpenPNM.Fluids.Water(network=pn)
water.regenerate()

#==============================================================================
'''Build Physics Objects'''
#==============================================================================
phys_water = OpenPNM.Physics.GenericPhysics(network=pn,fluid=water,name='standard_water_physics')
phys_water.add_method(prop='capillary_pressure',model='purcell',r_toroid=1e-5)
phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_water.add_method(prop='diffusive_conductance',model='bulk_diffusion')
phys_water.regenerate()

phys_air = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air,name='standard_air_physics')
phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_air.add_method(prop='diffusive_conductance',model='bulk_diffusion')
phys_air.regenerate()

#==============================================================================
'''Begin Simulations'''
#==============================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#------------------------------------------------------------------------------
#Initialize algorithm object
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(loglevel=20,loggername='OP',name='OP_1',network=pn)
a = pn.get_pore_indices(labels='bottom')
OP_1.run(invading_fluid='water',defending_fluid='air',inlets=a,npts=20)

#b = pn.get_pore_indices(labels='top')
#OP_1.evaluate_trapping(outlets=b)
#OP_1.plot_drainage_curve()

##-----------------------------------------------------------------------------
#'''Perform an Injection Experiment (InvasionPercolation)'''
##-----------------------------------------------------------------------------
##Initialize algorithm object
#IP_1 = OpenPNM.Algorithms.InvasionPercolation(loglevel=10,name='IP_1',network=pn)
#face = pn.get_pore_indices('right',indices=False)
#quarter = sp.rand(pn.num_pores(),)<.1
#inlets = pn.get_pore_indices()[face&quarter]
#outlets = pn.get_pore_indices('left')
#IP_1.run(invading_fluid=water,defending_fluid=air,inlets=inlets,outlets=outlets)

##----------------------------------------------------------------------
#'''Perform Fickian Diffusion'''
##----------------------------------------------------------------------
## Updating data based on the result of Percolation Algorithms
OP_1.update()
#IP_1.update()
##----------------------------------------------------------------------
## Initializing diffusion algorithm
Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg',network=pn)
## Assign Dirichlet boundary conditions
## BC1
BC1_pores = pn.get_pore_indices(labels='top')
Fickian_alg.set_pore_info(prop='Dirichlet',locations=BC1_pores,is_indices=True)
BC1_values = [0.8]
Fickian_alg.set_pore_data(subdomain='Dirichlet',prop='BCval',data=BC1_values,indices=BC1_pores)
## BC2
BC2_pores = pn.get_pore_indices(labels='bottom')
Fickian_alg.set_pore_info(prop='Dirichlet',locations=BC2_pores,is_indices=True)
BC2_values = [0.4]
Fickian_alg.set_pore_data(subdomain='Dirichlet',prop='BCval',data=BC2_values,indices=BC2_pores)
##----------------------------------------------------------------------
### Assign Neumann boundary conditions
### BC1
#BC1_pores = pn.get_pore_indices(subdomain='top')
#Fickian_alg.set_pore_info(prop='Dirichlet',locations=BC1_pores,is_indices=True)
#BC1_values = [0.5]
#Fickian_alg.set_pore_data(subdomain='Dirichlet',prop='BCval',data=BC1_values,indices=BC1_pores)
### BC2
#BC2_pores = pn.get_pore_indices(subdomain='bottom')
#Fickian_alg.set_pore_info(prop='Neumann_rate',locations=BC2_pores,is_indices=True)
#BC2_values = [2e-12]
#Fickian_alg.set_pore_data(subdomain='Neumann_rate',prop='BCval',data=BC2_values,indices=BC2_pores)
##----------------------------------------------------------------------
## Run simulation
Fickian_alg.run(active_fluid=air)
##----------------------------------------------------------------------
#
#
##Export to VTK
#OpenPNM.Visualization.VTK().write(pn,fluid=water)

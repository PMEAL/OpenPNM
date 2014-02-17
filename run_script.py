import OpenPNM

#==============================================================================
'''Build Topological Network'''
#==============================================================================
#pn = OpenPNM.Network.Cubic(name='cubic_1').generate(divisions=[15, 15, 15], lattice_spacing=[0.0001])
#pn = OpenPNM.Network.Delaunay(name='random_1',loglevel=10).generate(num_pores=1500,domain_size=[100,100,30])
#pn = OpenPNM.Network.Template(name='template_1',loglevel=10).generate(template=sp.ones((4,4,4),dtype=int),lattice_spacing=0.001)
pn = OpenPNM.Network.Sphere(name='sphere_1',loglevel=10).generate(radius=5,lattice_spacing=1)
#pn = OpenPNM.Network.TestNet()

#==============================================================================
'''Build Geometry'''
#==============================================================================
geom = OpenPNM.Geometry.Stick_and_Ball(network=pn, name='stick_and_ball', locations=pn.get_pore_indices())
geom.regenerate()

#==============================================================================
'''Build Fluids'''
#==============================================================================
air = OpenPNM.Fluids.Air(network=pn, init_cond={'temperature':300, 'pressure':100000})
air.apply_ICs(init_cond={'temperature':350, 'pressure':200000})  # experimental feature
air.regenerate()

water = OpenPNM.Fluids.Water(network=pn)
water.add_method(prop='diffusivity',prop_name='DAB',model='constant',value=1000000)
water.regenerate()

#==============================================================================
'''Build Physics Objects'''
#==============================================================================
phys_water = OpenPNM.Physics.GenericPhysics(network=pn, fluid=water, name='standard_water_physics')
phys_water.add_method(prop='capillary_pressure', model='purcell', r_toroid=1e-5)
phys_water.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_water.add_method(prop='diffusive_conductance', prop_name='gdAB', model='bulk_diffusion', diffusivity='DAB')
phys_water.regenerate()

phys_air = OpenPNM.Physics.GenericPhysics(network=pn, fluid=air, name='standard_air_physics')
phys_air.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_air.add_method(prop='diffusive_conductance', model='bulk_diffusion')
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
OP_1.update(Pc=3000)
###IP_1.update()
####----------------------------------------------------------------------
#### Initializing diffusion algorithm
Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loglevel=20, loggername='Fickian', name='Fickian_alg',network=pn)
## Assign Dirichlet boundary conditions
## BC1
BC1_pores = pn.get_pore_indices(labels='top')
Fickian_alg.set_pore_info(label='Dirichlet', locations=BC1_pores)
BC1_values = 0.8
Fickian_alg.set_pore_data(prop='BCval', data=BC1_values, locations=BC1_pores)
## BC2
BC2_pores = pn.get_pore_indices(labels='bottom')
Fickian_alg.set_pore_info(label='Dirichlet', locations=BC2_pores)
BC2_values = 0.4
Fickian_alg.set_pore_data(prop='BCval', data=BC2_values, locations=BC2_pores)
###----------------------------------------------------------------------
#### Assign Neumann boundary conditions
### BC1
#BC1_pores = pn.get_pore_indices(labels='top')
#Fickian_alg.set_pore_info(label='Dirichlet',locations=BC1_pores)
#BC1_values = 0.5
#Fickian_alg.set_pore_data(prop='BCval',data=BC1_values,locations=BC1_pores)
### BC2
#BC2_pores = pn.get_pore_indices(labels='bottom')
#Fickian_alg.set_pore_info(label='Neumann_rate',locations=BC2_pores)
#BC2_values = 2e-12
#Fickian_alg.set_pore_data(prop='BCval',data=BC2_values,locations=BC2_pores)
###----------------------------------------------------------------------
### Run simulation
Fickian_alg.run(active_fluid=air)
Fickian_alg.update()

#
#Export to VTK
OpenPNM.Visualization.VTK().write(net=pn, fluid=air)

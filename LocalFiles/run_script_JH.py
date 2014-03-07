import OpenPNM
import scipy as sp

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(name='cubic_1',loglevel=10).generate(divisions=[15, 15, 15], lattice_spacing=[0.0001])
#pn = OpenPNM.Network.Delaunay(name='random_1',loglevel=10).generate(num_pores=1500,domain_size=[100,100,30])
#pn = OpenPNM.Network.TestNet()

#==============================================================================
'''Build Geometry'''
#==============================================================================
ps = pn.get_pore_indices('boundary')
bndry = OpenPNM.Geometry.Boundary(network=pn, name='boundary', pnums=ps)
bndry.regenerate()

ps = pn.get_pore_indices('internal')
ts = pn.get_throat_indices('all')
geom = OpenPNM.Geometry.Stick_and_Ball(network=pn, name='domain', pnums=ps, tnums=ts)
geom.regenerate()

#==============================================================================
'''Build Fluids'''
#==============================================================================
air = OpenPNM.Fluids.Air(network=pn, loglevel=10,init_cond={'temperature':300, 'pressure':100000})
air.regenerate()

water = OpenPNM.Fluids.Water(network=pn,loglevel=10)
water.add_method(prop='diffusivity',prop_name='DAB',model='constant',value=5e-12)
water.regenerate()

#==============================================================================
'''Build Physics Objects'''
#==============================================================================
phys_water = OpenPNM.Physics.GenericPhysics(network=pn, fluid='water',name='phys_water')
phys_water.add_method(prop='capillary_pressure', model='purcell', r_toroid=1e-5)
phys_water.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_water.add_method(prop='diffusive_conductance', prop_name='gdAB', model='bulk_diffusion', diffusivity='DAB')
phys_water.regenerate()

phys_air = OpenPNM.Physics.GenericPhysics(network=pn, fluid=air, name='phys_air')
phys_air.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_air.add_method(prop='diffusive_conductance', model='bulk_diffusion')
phys_air.regenerate()

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
#OP_1.plot_drainage_curve()

#------------------------------------------------------------------------------
'''Perform Fickian Diffusion'''
#------------------------------------------------------------------------------
# Updating data based on the result of Percolation Algorithms
OP_1.update(Pc=8000)
# Initializing diffusion algorithm
Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loglevel=10, loggername='Fickian', name='Fickian_alg',network=pn)
# Assign Dirichlet boundary conditions to some of the surface pores
#BC1
BC1_pores = pn.get_pore_indices(labels=['top','boundary'],mode='intersection')
Fickian_alg.set_pore_info(label='Dirichlet', locations=BC1_pores)
BC1_values = 0.6
Fickian_alg.set_pore_data(prop='BCval', data=BC1_values, locations=BC1_pores)
# BC2
BC2_pores = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
Fickian_alg.set_pore_info(label='Dirichlet', locations=BC2_pores)
BC2_values = 0.2
Fickian_alg.set_pore_data(prop='BCval', data=BC2_values, locations=BC2_pores)
# Run simulation
Fickian_alg.run(active_fluid=air)
Fickian_alg.update()
#------------------------------------------------------------------------------
#Export to VTK
#OpenPNM.Visualization.VTK().write(net=pn, fluids=[air,water])  # Slow for now
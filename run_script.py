import OpenPNM
import scipy as sp

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(loglevel=20,name='net')
pn.generate(divisions=[20, 20, 20], lattice_spacing=[0.0001],add_boundaries=True)

#==============================================================================
'''Build Geometry'''
#==============================================================================
geom = OpenPNM.Geometry.Toray090(network=pn,name='geom')
geom.set_locations(pores=pn.pores('internal'),throats='all')

boun = pn.add_geometry(subclass='Boundary',name='boun')
boun.set_locations(pores=pn.pores('boundary'))

pn.regenerate_geometries()

#==============================================================================
'''Build Fluids'''
#==============================================================================
air = OpenPNM.Fluids.Air(network=pn, loglevel=20,name='air')
air.apply_conditions(temperature=350, pressure=200000)

water = OpenPNM.Fluids.Water(network=pn,loglevel=20,name='water')
water.add_property(prop='diffusivity',prop_name='DAB',model='constant',value=5e-12)

#Use Network's Fluid regeneration method
pn.regenerate_fluids()

#==============================================================================
'''Build Physics Objects'''
#==============================================================================
phys_water = OpenPNM.Physics.BasePhysics(network=pn, fluid=water,geometry=geom,name='physwater')

phys_air = OpenPNM.Physics.BasePhysics(network=pn, fluid=air,geometry=geom,name='physair')

#Use Network's Physics regeneration method
pn.regenerate_physics()

#==============================================================================
'''Begin Simulations'''
#==============================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#------------------------------------------------------------------------------
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(loglevel=20,network=pn)
a = pn.pores(labels=['bottom','boundary'],mode='intersection')
OP_1.setup(invading_fluid=water,defending_fluid=air,inlets=a,npts=20)
OP_1.run()
#OP_1.plot_drainage_curve()

#------------------------------------------------------------------------------
'''Perform Fickian Diffusion'''
#------------------------------------------------------------------------------
Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loglevel=20, network=pn)
# Assign Dirichlet boundary conditions to top and bottom surface pores
BC1_pores = pn.pores(labels=['top','front'],mode='intersection')
Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)

BC2_pores = pn.pores(labels=['top','back'],mode='intersection')
Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=BC2_pores)

# Updating data based on the result of Percolation Algorithms
OP_1.update(Pc=11000)
# Run simulation
Fickian_alg.run(active_fluid=air)
Fickian_alg.update()

#------------------------------------------------------------------------------
'''Export to VTK'''
#------------------------------------------------------------------------------
#OpenPNM.Visualization.Vtp.write(filename='test.vtp',fluids=[air,water],network=pn)
vis = OpenPNM.Visualization.VTK()
vis.write(network=pn,fluids=[air,water])







import OpenPNM
import scipy as sp
import OpenPNM.Utilities.misc as misc

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(loglevel=30,name='net')
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
air = OpenPNM.Fluids.Air(network=pn, loglevel=30,name='air')
air.apply_conditions(temperature=350, pressure=200000)
air.add_property(prop='thermal_conductivity',model='constant',value=0.0262)
air.add_property(prop='electrical_conductivity',model='constant',value=1)

water = OpenPNM.Fluids.Water(network=pn,loglevel=30,name='water')
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
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(loglevel=30,network=pn)
a = pn.pores(labels=['bottom','boundary'],mode='intersection')
OP_1.setup(invading_fluid=water,defending_fluid=air,inlets=a,npts=100)
OP_1.run()

# Updating data based on the result of Percolation Algorithms
OP_1.update(Pc=11000)
#OP_1.plot_drainage_curve()
inlets = pn.get_pore_indices(labels = ['bottom'])
outlets = pn.get_pore_indices(labels = ['top'])


IP_1 = OpenPNM.Algorithms.InvasionPercolation(network = pn, name = 'OP_1', loglevel = 30)
IP_1.setup(invading_fluid = water, defending_fluid = air, inlets = inlets, outlets = outlets, end_condition = 'breakthrough')
IP_1.run()
IP_1.update()
#------------------------------------------------------------------------------
'''Perform Fickian Diffusion'''
#------------------------------------------------------------------------------
alg = OpenPNM.Algorithms.FickianDiffusion(loglevel=10, network=pn)
# Assign Dirichlet boundary conditions to top and bottom surface pores
BC1_pores = pn.pores(labels=['top','boundary'],mode='intersection')
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=sp.log(1-.6), pores=BC1_pores)

BC2_pores = pn.pores(labels=['bottom','boundary'],mode='intersection')
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=sp.log(1-0.4), pores=BC2_pores)

phys_air.add_property(prop='multiphase',model='conduit_conductance',
                      conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='loose')
pn.regenerate_physics()

# Run simulation
alg.setup(conductance = 'conduit_diffusive_conductance',fluid=air)
alg.run()
alg.update()
Deff = alg.calc_eff_diffusivity()
print(Deff)

#Deff = OpenPNM.Algorithms.EffectiveProperty(network=pn)
#Deff.setup(algorithm=alg,fluid=air,conductance='conduit_diffusive_conductance',quantity='mole_fraction')
#a = Deff.run(algorithm = alg)
#


#------------------------------------------------------------------------------
'''Export to VTK'''
#------------------------------------------------------------------------------

#OpenPNM.Visualization.Vtp.write(filename='test.vtp',fluids=[air,water],network=pn)

vis = OpenPNM.Visualization.VTK()
vis.write(network=pn,fluids=[air,water])

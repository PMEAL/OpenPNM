import OpenPNM
import scipy as sp

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(loglevel=30,name='net')
pn.generate(divisions=[5, 5, 5], lattice_spacing=[0.0001],add_boundaries=True)

#==============================================================================
'''Build Geometry'''
#==============================================================================
#Example of assigning pores/throats AFTER instantiation
geom = OpenPNM.Geometry.Toray090(network=pn)
Ps = pn.pores('boundary',mode='difference')
Ts = pn.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
geom.set_locations(pores=Ps,throats=Ts)
geom.generate()

#Example of assigning pores/throats DURING initialization.
Ps = pn.pores('boundary')
Ts = pn.find_neighbor_throats(pores=Ps,mode='not_intersection')
boun = OpenPNM.Geometry.Boundary(network=pn,pores=Ps,throats=Ts)
boun.generate()

#==============================================================================
'''Build Fluids'''
#==============================================================================
#Fluids exist everywhere so don't need to be given pores/throats
air = OpenPNM.Fluids.Air(network=pn)
air['pore.thermal_conductivity'] = 0.0262
air['pore.electrical_conductivity'] = 1
air.generate()

water = OpenPNM.Fluids.Water(network=pn)
water['pore.DAB'] = 5e-12
water.generate()

#==============================================================================
'''Build Physics'''
#==============================================================================
#Example of assigning pores/throats DURING initialization.
phys_water = OpenPNM.Physics.Standard(network=pn,fluid=water,pores=pn.pores(),throats=pn.throats())
phys_water.generate()

#Example of assigning pores/throats AFTER instantiation
phys_air = OpenPNM.Physics.Standard(network=pn,fluid=air)
phys_air.set_locations(pores=pn.pores(),throats=pn.throats())
phys_air.generate()

#==============================================================================
'''Begin Simulations'''
#==============================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#------------------------------------------------------------------------------
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn)
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
                      conductance = 'diffusive_conductance', prop_name='conduit_diffusive_conductance',mode='strict')
pn.regenerate_physics()

# Run simulation
alg.setup(conductance = 'diffusive_conductance',fluid=air)
alg.run()
alg.update()
Deff = alg.calc_eff_diffusivity()


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

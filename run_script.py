import OpenPNM

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(name='net')
pn.generate(divisions=[5, 5, 5], lattice_spacing=[0.0001],add_boundaries=True)

#==============================================================================
'''Build Geometry'''
#==============================================================================
Ps = pn.pores('boundary',mode='difference')
Ts = pn.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
geom = OpenPNM.Geometry.Toray090(network=pn,pores=Ps,throats=Ts,dynamic_data=True)

Ps = pn.pores('boundary')
Ts = pn.find_neighbor_throats(pores=Ps,mode='not_intersection')
boun = OpenPNM.Geometry.Boundary(network=pn,pores=Ps,throats=Ts)

#==============================================================================
'''Build Fluids'''
#==============================================================================
#Fluids exist everywhere so don't need to be given pores/throats
air = OpenPNM.Fluids.Air(network=pn)
air['pore.Dac'] = 1e-7  # Add custom properties directly
water = OpenPNM.Fluids.Water(network=pn)

#==============================================================================
'''Build Physics'''
#==============================================================================
Ps = pn.pores()
Ts = pn.throats()
phys_water = OpenPNM.Physics.Standard(network=pn,fluid=water,pores=Ps,throats=Ts)
phys_air = OpenPNM.Physics.Standard(network=pn,fluid=air,pores=Ps,throats=Ts)
#Add some additional models to phys_air
phys_air.add_model(model=OpenPNM.Physics.models.diffusive_conductance.bulk_diffusion,
                   propname='throat.gdiff_ac',
                   pore_diffusivity='pore.Dac')

#==============================================================================
'''Begin Simulations'''
#==============================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#------------------------------------------------------------------------------
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn,loglevel=20)
Ps = pn.pores(labels=['bottom_face'])
OP_1.setup(invading_fluid=water,defending_fluid=air,inlets=Ps)
OP_1.run()
OP_1.update(Pc=7000)
#OP_1.plot_drainage_curve()

#------------------------------------------------------------------------------
'''Perform Invasion Percolation'''
#------------------------------------------------------------------------------
inlets = pn.get_pore_indices(labels = ['bottom_face'])
outlets = pn.get_pore_indices(labels = ['top_face'])
IP_1 = OpenPNM.Algorithms.InvasionPercolation(network = pn, name = 'IP_1', loglevel = 30)
IP_1.setup(invading_fluid = water, defending_fluid = air, inlets = inlets, outlets = outlets, end_condition = 'breakthrough')
IP_1.run()
IP_1.update()

#------------------------------------------------------------------------------
'''Perform Fickian Diffusion'''
#------------------------------------------------------------------------------
alg = OpenPNM.Algorithms.FickianDiffusion(loglevel=10, network=pn)
# Assign Dirichlet boundary conditions to top and bottom surface pores
BC1_pores = pn.pores(labels=['top_face'])
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
BC2_pores = pn.pores(labels=['bottom_face'])
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
#Add new model to air's physics that accounts for water occupancy
phys_air.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
                   propname='throat.conduit_diffusive_conductance',
                   throat_conductance='throat.diffusive_conductance')
#Use newly defined diffusive_conductance in the diffusion calculation
alg.setup(conductance = 'throat.conduit_diffusive_conductance',fluid=air)
alg.run()
alg.update()
Deff = alg.calc_eff_diffusivity()

#------------------------------------------------------------------------------
'''Export to VTK'''
#------------------------------------------------------------------------------
vis = OpenPNM.Visualization.VTK()
vis.write(network=pn,fluids=[air,water])

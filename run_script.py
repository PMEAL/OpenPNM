import OpenPNM as op
print('-----> Using OpenPNM version: '+op.__version__)

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = op.Network.Cubic(shape=[5,6,7],spacing=0.0001,name='net')
pn.add_boundaries()

#==============================================================================
'''Build Geometry'''
#==============================================================================
Ps = pn.pores('boundary',mode='not')
Ts = pn.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
geom = op.Geometry.Toray090(network=pn,pores=Ps,throats=Ts)

Ps = pn.pores('boundary')
Ts = pn.find_neighbor_throats(pores=Ps,mode='not_intersection')
boun = op.Geometry.Boundary(network=pn,pores=Ps,throats=Ts)

#==============================================================================
'''Build Phases'''
#==============================================================================
air = op.Phases.Air(network=pn,name='air')
air['pore.Dac'] = 1e-7  # Add custom properties directly
water = op.Phases.Water(network=pn,name='water')

#==============================================================================
'''Build Physics'''
#==============================================================================
Ps = pn.pores()
Ts = pn.throats()
phys_water = op.Physics.Standard(network=pn,phase=water,pores=Ps,throats=Ts)
phys_air = op.Physics.Standard(network=pn,phase=air,pores=Ps,throats=Ts)
#Add some additional models to phys_air
phys_air.models.add(model=op.Physics.models.diffusive_conductance.bulk_diffusion,
                    propname='throat.gdiff_ac',
                    pore_diffusivity='pore.Dac')

#==============================================================================
'''Begin Simulations'''
#==============================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#------------------------------------------------------------------------------
OP_1 = op.Algorithms.OrdinaryPercolation(network=pn,invading_phase=water,defending_phase=air)
Ps = pn.pores(labels=['bottom_boundary'])
OP_1.run(inlets=Ps)
OP_1.return_results(Pc=7000)

#------------------------------------------------------------------------------
'''Perform Invasion Percolation'''
#------------------------------------------------------------------------------
inlets = pn.pores('bottom_boundary')
outlets = pn.pores('top_boundary')
IP_1 = op.Algorithms.InvasionPercolation(network=pn,name='IP_1')
IP_1.run(phase=water,inlets=inlets)
IP_1.return_results()

#------------------------------------------------------------------------------
'''Perform Fickian Diffusion'''
#------------------------------------------------------------------------------
alg = op.Algorithms.FickianDiffusion(network=pn,phase=air)
# Assign Dirichlet boundary conditions to top and bottom surface pores
BC1_pores = pn.pores('right_boundary')
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
BC2_pores = pn.pores('left_boundary')
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
#Add new model to air's physics that accounts for water occupancy
phys_air.models.add(model=op.Physics.models.multiphase.conduit_conductance,
                    propname='throat.conduit_diffusive_conductance',
                    throat_conductance='throat.diffusive_conductance',
                    throat_occupancy='throat.occupancy',
                    pore_occupancy='pore.occupancy',
                    mode='strict',
                    factor=0)
#Use desired diffusive_conductance in the diffusion calculation (conductance for the dry network or water-filled network)
alg.run(conductance='throat.diffusive_conductance')
alg.return_results()
Deff = alg.calc_eff_diffusivity()

try:
    # this creates a time step x num_pores, which is what the animated object needs
    inv_seq = water['pore.IP_inv_seq'].squeeze()
    history = []
    for i in sorted(set(inv_seq)):
      history.append( (inv_seq != 0) & (inv_seq < i) )

except Exception as e:
    pass

#------------------------------------------------------------------------------
'''Export to VTK'''
#------------------------------------------------------------------------------
op.export()

import OpenPNM
print('-----> Using OpenPNM version: '+OpenPNM.__version__)

ctrl = OpenPNM.Base.Controller()
#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(shape=[5,6,7],spacing=0.0001,name='net')
pn.add_boundaries()

#==============================================================================
'''Build Geometry'''
#==============================================================================
Ps = pn.pores('boundary',mode='not')
Ts = pn.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
geom = OpenPNM.Geometry.Toray090(network=pn,pores=Ps,throats=Ts)
geom.models['pore.seed']['seed'] = 0
geom.models['pore.seed']['regen_mode'] = 'normal'
geom.regenerate()

Ps = pn.pores('boundary')
Ts = pn.find_neighbor_throats(pores=Ps,mode='not_intersection')
boun = OpenPNM.Geometry.Boundary(network=pn,pores=Ps,throats=Ts)

#==============================================================================
'''Build Phases'''
#==============================================================================
air = OpenPNM.Phases.Air(network=pn,name='air')
water = OpenPNM.Phases.Water(network=pn,name='water')

#==============================================================================
'''Build Physics'''
#==============================================================================
Ps = pn.pores()
Ts = pn.throats()
phys_water = OpenPNM.Physics.Standard(network=pn,phase=water,pores=Ps,throats=Ts)
phys_air = OpenPNM.Physics.Standard(network=pn,phase=air,pores=Ps,throats=Ts)
#Add some additional models to phys_air
phys_air.models.add(model=OpenPNM.Physics.models.capillary_pressure.static_pressure,
                    propname='pore.static_pressure',
                    regen_mode='deferred')

#==============================================================================
'''Begin Simulations'''
#==============================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#------------------------------------------------------------------------------
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn,invading_phase=water,defending_phase=air)
Ps = pn.pores(labels=['bottom_boundary'])
OP_1.run(inlets=Ps)
OP_1.return_results(Pc=7000)

#------------------------------------------------------------------------------
'''Perform Invasion Percolation'''
#------------------------------------------------------------------------------
IP_1 = OpenPNM.Algorithms.InvasionPercolation(network=pn)
inlets = pn.pores('front_boundary')
outlets = pn.pores('back_boundary')
IP_1.run(phase=water,inlets=inlets)
IP_1.apply_flow(flowrate=1e-15)
IP_1.return_results()

#------------------------------------------------------------------------------
'''Perform Invasion Percolation using Version 2'''
#------------------------------------------------------------------------------
IP_2 = OpenPNM.Algorithms.InvasionPercolation2(network=pn)
inlets = pn.pores('front_boundary')
IP_2.setup(phase=water, inlets=inlets)
filled = False
while not filled:
    filled = IP_2.run(nsteps=10)

#------------------------------------------------------------------------------
'''Perform Fickian Diffusion'''
#------------------------------------------------------------------------------
alg = OpenPNM.Algorithms.FickianDiffusion(network=pn,phase=air)
# Assign Dirichlet boundary conditions to top and bottom surface pores
BC1_pores = pn.pores('right_boundary')
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
BC2_pores = pn.pores('left_boundary')
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
#Add new model to air's physics that accounts for water occupancy
phys_air.models.add(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
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

#------------------------------------------------------------------------------
'''Export to VTK'''
#------------------------------------------------------------------------------
ctrl.export()


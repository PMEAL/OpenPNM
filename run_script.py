import OpenPNM

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(name='net',loglevel=20,divisions=[5,5,5],lattice_spacing=[0.0001],add_boundaries=True)

#==============================================================================
'''Build Geometry'''
#==============================================================================
Ps = pn.pores('boundary',mode='difference')
Ts = pn.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
geom = OpenPNM.Geometry.Toray090(network=pn,pores=Ps,throats=Ts)

Ps = pn.pores('boundary')
Ts = pn.find_neighbor_throats(pores=Ps,mode='not_intersection')
boun = OpenPNM.Geometry.Boundary(network=pn,pores=Ps,throats=Ts)

#==============================================================================
'''Build Fluids'''
#==============================================================================
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
OP_1.run(invading_fluid=water,defending_fluid=air,inlets=Ps)
OP_1.update(Pc=7000)
#OP_1.plot_drainage_curve()

#------------------------------------------------------------------------------
'''Perform a Drainage Experiment on a SUB-network'''
#------------------------------------------------------------------------------
#Create a sub-network
#import OpenPNM.Utilities.Subsets as subs
#sub_pn = subs.subset_network(pn,pores=pn.pores(geom.name))
#sub_water = subs.subset_fluid(fluid=water,subnet=sub_pn)
##Run standard algorithm on subnet, and subfluid
#OP_2 = OpenPNM.Algorithms.OrdinaryPercolation(network=sub_pn)
#OP_2.setup(invading_fluid=sub_water,inlets=sub_pn.pores('bottom'))
#OP_2.run()

#------------------------------------------------------------------------------
'''Perform Invasion Percolation'''
##------------------------------------------------------------------------------
inlets = pn.pores('bottom_face')
outlets = pn.pores('top_face')
IP_1 = OpenPNM.Algorithms.InvasionPercolation(network = pn, name = 'IP_1', loglevel = 30)
IP_1.run(invading_fluid = water, defending_fluid = air, inlets = inlets, outlets = outlets, end_condition = 'breakthrough')
IP_1.update()

#------------------------------------------------------------------------------
'''Perform Fickian Diffusion'''
#------------------------------------------------------------------------------
alg = OpenPNM.Algorithms.FickianDiffusion(loglevel=20, network=pn)
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
alg.run(conductance='throat.diffusive_conductance',fluid=air)
alg.update()
Deff = alg.calc_eff_diffusivity()

# this creates a time step x num_pores, which is what the animated object needs
inv_seq = water['pore.IP_inv_seq'].squeeze()
history = []
for i in sorted(set(inv_seq)):
  history.append( (inv_seq != 0) & (inv_seq < i) )

try:
  # try to perofrm an animated 3D rendering
  from OpenPNM.Postprocessing.Graphics import Scene, Wires
  wires = Wires(pn['pore.coords'], pn['throat.conns'], history)
  scene = Scene()
  scene.add_actors([wires])
  scene.play()

except Exception as e:
  #------------------------------------------------------------------------------
  '''Export to VTK'''
  #------------------------------------------------------------------------------
  import OpenPNM.Postprocessing.VTK as vtk
  vtk.save(network=pn,fluids=[air,water])

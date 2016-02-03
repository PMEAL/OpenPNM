import OpenPNM as op
import time
st = time.time()

# =============================================================================
'''Build Topological Network'''
# =============================================================================
pn = op.Network.Cubic(shape=[5, 6, 7], spacing=0.0001, name='net')
pn.add_boundaries()

# =============================================================================
'''Build Geometry'''
# =============================================================================
Ps = pn.pores('boundary', mode='not')
Ts = pn.find_neighbor_throats(pores=Ps, mode='intersection', flatten=True)
geom = op.Geometry.Toray090(network=pn, pores=Ps, throats=Ts)
Ps = pn.pores('boundary')
Ts = pn.find_neighbor_throats(pores=Ps, mode='not_intersection')
boun = op.Geometry.Boundary(network=pn, pores=Ps, throats=Ts)

# =============================================================================
'''Build Phases'''
# =============================================================================
air = op.Phases.Air(network=pn)
air['pore.Dac'] = 1e-7  # Add custom properties directly
water = op.Phases.Water(network=pn)

# =============================================================================
'''Build Physics'''
# =============================================================================
Ps = pn.pores()
Ts = pn.throats()
phys_water = op.Physics.Standard(network=pn, phase=water, pores=Ps, throats=Ts)
phys_air = op.Physics.Standard(network=pn, phase=air, pores=Ps, throats=Ts)

# =============================================================================
'''Begin Simulations'''
# =============================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
# -----------------------------------------------------------------------------
drainage = op.Algorithms.Drainage(network=pn)
drainage.setup(invading_phase=water)
drainage.set_inlets(pores=pn.pores(labels=['bottom_boundary']))
drainage.run()
drainage.return_results(Pc=7000)

# -----------------------------------------------------------------------------
'''Perform Fickian Diffusion in Partially Invaded Network'''
# -----------------------------------------------------------------------------
alg = op.Algorithms.FickianDiffusion(network=pn, phase=air)
# Assign Dirichlet boundary conditions to top and bottom surface pores
BC1_pores = pn.pores('right_boundary')
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
BC2_pores = pn.pores('left_boundary')
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
# Add new model to air's physics that accounts for water occupancy
phys_air.models.add(model=op.Physics.models.multiphase.conduit_conductance,
                    propname='throat.conduit_diffusive_conductance',
                    throat_conductance='throat.diffusive_conductance',
                    throat_occupancy='throat.occupancy',
                    pore_occupancy='pore.occupancy',
                    mode='strict',
                    factor=0)
# Use desired conduit_diffusive_conductance in the diffusion calculation
alg.run(conductance='throat.conduit_diffusive_conductance')
alg.return_results()
Deff = alg.calc_eff_diffusivity()

# -----------------------------------------------------------------------------
'''Export to VTK'''
# -----------------------------------------------------------------------------
op.export()
print("sim time:" + str(time.time()-st))

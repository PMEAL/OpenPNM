import OpenPNM as op

# Build Topological Network
pn = op.Network.Cubic(shape=[15, 15, 15], spacing=1, name='net')
pn.add_boundaries()

# Build Geometry objects
Ps = pn.pores('boundary', mode='not')
Ts = pn.find_neighbor_throats(pores=Ps, mode='intersection', flatten=True)
geom = op.Geometry.Stick_and_Ball(network=pn, pores=Ps, throats=Ts)
Ps = pn.pores('boundary')
Ts = pn.find_neighbor_throats(pores=Ps, mode='not_intersection')
boun = op.Geometry.Boundary(network=pn, pores=Ps, throats=Ts)

# Build Phases
air = op.Phases.Air(network=pn, name='air')
water = op.Phases.Water(network=pn, name='water')

# Build Physics
Ps = pn.pores()
Ts = pn.throats()
phys_water = op.Physics.Standard(network=pn, phase=water, pores=Ps, throats=Ts)
phys_air = op.Physics.Standard(network=pn, phase=air, pores=Ps, throats=Ts)

# Perform a Drainage Experiment (OrdinaryPercolation)
MIP = op.Algorithms.Drainage(network=pn)
MIP.setup(invading_phase=water, defending_phase=air)
Ps = pn.pores(labels='bottom_boundary')
MIP.set_inlets(pores=Ps)
MIP.run()
MIP.return_results(Pc=7000)

# Perform Invasion Percolation
IP = op.Algorithms.InvasionPercolation(network=pn)
IP.setup(phase=water)
Ps = pn.pores('bottom_boundary')
IP.set_inlets(pores=Ps)
IP.run()
IP.return_results()

# Perform Fickian Diffusion
alg = op.Algorithms.FickianDiffusion(network=pn, phase=air)
BC1_pores = pn.pores('right_boundary')
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
BC2_pores = pn.pores('left_boundary')
alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
alg.run()
alg.return_results()
Deff = alg.calc_eff_diffusivity()
print(Deff/air['pore.diffusivity'][0])

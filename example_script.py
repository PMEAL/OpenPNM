import OpenPNM as op

# Create Topological Network object
pn = op.Network.Cubic(shape=[15, 5, 5], spacing=1)
pn.add_boundaries()

# Create Geometry object for internal pores
Ps = pn.pores('boundary', mode='not')
Ts = pn.find_neighbor_throats(pores=Ps, mode='intersection', flatten=True)
geom = op.Geometry.Stick_and_Ball(network=pn, pores=Ps, throats=Ts)
# Create Geometry object for boundary pores
Ps = pn.pores('boundary')
Ts = pn.find_neighbor_throats(pores=Ps, mode='not_intersection')
boun = op.Geometry.Boundary(network=pn, pores=Ps, throats=Ts)

# Create Phase objects
air = op.Phases.Air(network=pn, name='air')
water = op.Phases.Water(network=pn, name='water')

# Create one Physics object for each phase
Ps = pn.pores()
Ts = pn.throats()
phys_water = op.Physics.Standard(network=pn, phase=water, pores=Ps, throats=Ts)
phys_air = op.Physics.Standard(network=pn, phase=air, pores=Ps, throats=Ts)

# Begin Simulations
# 1. Perform a Drainage Experiment
drainage = op.Algorithms.Drainage(network=pn)
drainage.setup(invading_phase=water, defending_phase=air)
Ps = pn.pores(labels=['bottom_boundary'])
drainage.set_inlets(pores=Ps)
drainage.run()
drainage.return_results(Pc=7000)

# 2. Perform Invasion Percolation
inlets = pn.pores('bottom_boundary')
IP_1 = op.Algorithms.InvasionPercolation(network=pn)
IP_1.run(phase=water, inlets=inlets)
IP_1.return_results()

# 3. Perform Fickian Diffusion'''
fickian = op.Algorithms.FickianDiffusion(network=pn, phase=air)
Ps_bc1 = pn.pores('right_boundary')
fickian.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=Ps_bc1)
Ps_bc2 = pn.pores('left_boundary')
fickian.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=Ps_bc2)
fickian.run()
fickian.return_results()

# Export to VTK
op.export_data(network=pn, filename='example', fileformat='VTK')

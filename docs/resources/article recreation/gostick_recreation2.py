import OpenPNM

Lc = 40.5e-6 #Lattice constant used in [1] for SGL 10BA
#set up network "sgl"
sgl = OpenPNM.Network.Cubic([26, 26, 10], spacing=Lc, name='sgl')
sgl.add_boundaries()

#set up geometries, "geo" and "boun"
Ps = sgl.pores('boundary',mode='difference')
Ts = sgl.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
geo = OpenPNM.Geometry.SGL10(network=sgl,pores=Ps,throats=Ts,name='geo')

Ps = sgl.pores('boundary')
Ts = sgl.find_neighbor_throats(pores=Ps,mode='not_intersection')
boun = OpenPNM.Geometry.Boundary(network=sgl,pores=Ps,throats=Ts,name='boun')

#constrict throats in the y and z directions
throats = sgl.throats('geo')
connected_pores = sgl.find_connected_pores(throats)
x1 = [sgl['pore.coords'][pair[0]][0] for pair in connected_pores]
x2 = [sgl['pore.coords'][pair[1]][0] for pair in connected_pores]
same_x = [x - y == 0 for x, y in zip(x1,x2)]
factor = [s*.95 + (not s)*1 for s in same_x]
throat_diameters = sgl['throat.diameter'][throats]*factor

#remove the regeneration ability of the diameter pore and throat properties
geo.remove_model('pore.diameter')
geo.remove_model('throat.diameter')
boun.remove_model('pore.diameter')
boun.remove_model('throat.diameter')

#reset aspects relying on pore and throat sizes
geo.regenerate()
boun.regenerate()

import OpenPNM.Utilities.IO as io
io.VTK.save(network=sgl)

#set up phases
air = OpenPNM.Phases.Air(network = sgl, name = 'air')
water = OpenPNM.Phases.Water(network = sgl, name = 'water')

#reset pore contact angle
water['pore.contact_angle'] = 100
#remove the
water.remove_model('pore.contact_angle')

#create physics objects associated with our phases
Ps = sgl.pores()
Ts = sgl.throats()
phys_water = OpenPNM.Physics.Standard(network=sgl,phase=water,pores=Ps,throats=Ts,dynamic_data=True,name='standard_water_physics')
phys_air = OpenPNM.Physics.Standard(network=sgl,phase=air,pores=Ps,throats=Ts,dynamic_data=True,name='standard_air_physics')

inlets = sgl.pores('bottom_boundary')
used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]

#using every other pore in the bottom and boundary as an inlet
#prevents extremely small diffusivity and permeability values in the z direction
used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]

OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=sgl,invading_phase = water, defending_phase = air, )
OP_1.run(inlets = used_inlets,npts=100)

#Update the simulation until saturation is at 50%
OP_1.return_results(sat=0.5)

#adding multiphase conductances
phys_air.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
           propname='throat.conduit_diffusive_conductance',
           throat_conductance='throat.diffusive_conductance')
phys_water.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
           propname='throat.conduit_diffusive_conductance',
           throat_conductance='throat.diffusive_conductance')
phys_air.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
           propname='throat.conduit_hydraulic_conductance',
           throat_conductance='throat.hydraulic_conductance')
phys_water.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
           propname='throat.conduit_hydraulic_conductance',
           throat_conductance='throat.hydraulic_conductance')
           
#setting up the 8 StokesFlow and FickianDiffusion algorithms
Stokes_alg_single_phase_air = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_single_phase_air',network=sgl,phase=air)
Stokes_alg_single_phase_water = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_single_phase_water',network=sgl,phase=water)

Fickian_alg_single_phase_air = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_single_phase_air',network=sgl,phase=air)
Fickian_alg_single_phase_water = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_single_phase_water',network=sgl,phase=water)

Stokes_alg_multi_phase_air = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_multi_phase_air',network=sgl,phase=air)
Stokes_alg_multi_phase_water = OpenPNM.Algorithms.StokesFlow(name='Stokes_alg_multi_phase_water',network=sgl,phase=water)

Fickian_alg_multi_phase_air = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_multi_phase_air',network=sgl,phase=air)
Fickian_alg_multi_phase_water = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg_multi_phase_water',network=sgl,phase=water)

#setting boundary conditions
BC1_pores = sgl.pores(labels='bottom_boundary')
BC2_pores = sgl.pores(labels='top_boundary')

#BC1
Stokes_alg_single_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.6,pores=BC1_pores)
Stokes_alg_single_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.6,pores=BC1_pores)
Fickian_alg_single_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=.6,pores=BC1_pores)
Fickian_alg_single_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=.6,pores=BC1_pores)

Stokes_alg_multi_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.6,pores=BC1_pores)
Stokes_alg_multi_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.6,pores=BC1_pores)
Fickian_alg_multi_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=.6,pores=BC1_pores)
Fickian_alg_multi_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=.6,pores=BC1_pores)

#BC2
Stokes_alg_single_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.2,pores=BC2_pores)
Stokes_alg_single_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.2,pores=BC2_pores)
Fickian_alg_single_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=.2,pores=BC2_pores)
Fickian_alg_single_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=.2,pores=BC2_pores)

Stokes_alg_multi_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.2,pores=BC2_pores)
Stokes_alg_multi_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.2,pores=BC2_pores)
Fickian_alg_multi_phase_air.set_boundary_conditions(bctype='Dirichlet',bcvalue=.2,pores=BC2_pores)
Fickian_alg_multi_phase_water.set_boundary_conditions(bctype='Dirichlet',bcvalue=.2,pores=BC2_pores)

#run algorithms with proper conduit conductance
Stokes_alg_single_phase_air.run(conductance = 'hydraulic_conductance')
Stokes_alg_single_phase_water.run(conductance = 'hydraulic_conductance')
Fickian_alg_single_phase_air.run(conductance = 'diffusive_conductance')
Fickian_alg_single_phase_water.run(conductance = 'diffusive_conductance')

Stokes_alg_multi_phase_air.run(conductance = 'conduit_hydraulic_conductance')
Stokes_alg_multi_phase_water.run(conductance = 'conduit_hydraulic_conductance')
Fickian_alg_multi_phase_air.run(conductance = 'conduit_diffusive_conductance')
Fickian_alg_multi_phase_water.run(conductance = 'conduit_diffusive_conductance')

#calc effective properties
effective_permeability_air_single = Stokes_alg_single_phase_air.calc_eff_permeability()
effective_diffusivity_air_single = Fickian_alg_single_phase_air.calc_eff_diffusivity()
effective_permeability_water_single = Stokes_alg_single_phase_water.calc_eff_permeability()
effective_diffusivity_water_single = Fickian_alg_single_phase_water.calc_eff_diffusivity()

effective_permeability_air_multi = Stokes_alg_multi_phase_air.calc_eff_permeability()
effective_diffusivity_air_multi = Fickian_alg_multi_phase_air.calc_eff_diffusivity()
effective_permeability_water_multi = Stokes_alg_multi_phase_water.calc_eff_permeability()
effective_diffusivity_water_multi = Fickian_alg_multi_phase_water.calc_eff_diffusivity()

relative_eff_perm_air = effective_permeability_air_multi/effective_permeability_air_single
relative_eff_perm_water = effective_permeability_water_multi/effective_permeability_water_single
relative_eff_diff_air = effective_diffusivity_air_multi/effective_diffusivity_air_single
relative_eff_diff_water = effective_diffusivity_water_multi/effective_diffusivity_water_single


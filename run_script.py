import OpenPNM
import scipy as sp

#======================================================================
'''Build Topological Network'''
#======================================================================
pn = OpenPNM.Network.Cubic(name='cubic_1').generate(divisions=[35,35,35],lattice_spacing=[0.0001])

#======================================================================
'''Build Geometry'''
#======================================================================
geom = OpenPNM.Geometry.GenericGeometry(loglevel=10,network=pn,name='stick_and_ball')
geom.add_method(prop='pore_seed',model='random')
geom.add_method(prop='throat_seed',model='neighbor_min')
geom.add_method(prop='pore_diameter',model='sphere',name='weibull_min',shape=2.5,loc=6e-6,scale=2e-5)
geom.add_method(prop='throat_diameter',model='cylinder',name='weibull_min',shape=2.5,loc=6e-6,scale=2e-5)
geom.add_method(prop='pore_volume',model='sphere')
geom.add_method(prop='throat_volume',model='cylinder')
geom.add_method(prop='throat_length',model='straight')
geom.regenerate()

#======================================================================
'''Build Fluids'''
#======================================================================
air = OpenPNM.Fluids.GenericFluid(loggername='AIR',loglevel=10,network=pn,name='air')
air.set_pore_data(prop='Pc',data=132.65)
air.set_pore_data(prop='Tc',data=3.771e6)
air.set_pore_data(prop='MW',data=0.0291)
air.add_method(prop='diffusivity',model='Fuller',MA=0.03199,MB=0.0291,vA=16.3,vB=19.7)
air.add_method(prop='viscosity',model='Reynolds',uo=0.001,b=0.1)
air.add_method(prop='molar_density',model='ideal_gas',R=8.314)
air.regenerate()

water = OpenPNM.Fluids.GenericFluid(loggername='AIR',loglevel=10,network=pn,name='water')
water.set_pore_data(prop='Pc',data=132.65)
water.set_pore_data(prop='Tc',data=3.771e6)
water.set_pore_data(prop='MW',data=0.0291)
water.add_method(prop='diffusivity',model='constant',value=1e-12)
water.add_method(prop='viscosity',model='constant',value=0.001)
water.add_method(prop='molar_density',model='constant',value=44445)
water.add_method(prop='surface_tension',model='constant',value=0.072)
water.add_method(prop='contact_angle',model='constant',value=110)
water.regenerate()

#======================================================================
'''Build Physics Objects'''
#======================================================================
phys_water = OpenPNM.Physics.GenericPhysics(network=pn,fluid=water,name='standard_water_physics')
phys_water.add_method(prop='capillary_pressure',model='purcell',r_toroid=1e-5)
phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_water.add_method(prop='diffusive_conductance',model='bulk_diffusion')
phys_water.regenerate()
phys_air = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air,name='standard_air_physics')
phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_air.add_method(prop='diffusive_conductance',model='bulk_diffusion')
phys_air.regenerate()

#======================================================================
'''Begin Simulations'''
#======================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#----------------------------------------------------------------------
#Initialize algorithm object
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(loglevel=10,loggername="OP",name='OP_1',network=pn)
a = pn.get_pore_indices(subdomain='bottom')
#Run algorithm
OP_1.run(invading_fluid='water',defending_fluid='air',inlets=a,npts=20)

#b = pn.get_pore_indices(subdomain='top')
#OP_1.evaluate_trapping(outlets=b)
#OP_1.plot_drainage_curve()

##----------------------------------------------------------------------
#'''Perform an Injection Experiment (InvasionPercolation)'''
##----------------------------------------------------------------------
##Initialize algorithm object
#IP_1 = OpenPNM.Algorithms.InvasionPercolation()
##Apply desired/necessary pore scale physics methods
#OpenPNM.Physics.CapillaryPressure.Washburn(pn,water2)
#face = pn.pore_properties['type']==3
#quarter = sp.rand(pn.get_num_pores(),)<.1
#inlets = pn.pore_properties['numbering'][face&quarter]
#outlets = pn.pore_properties['numbering'][pn.pore_properties['type']==4]
#IP_1.run(pn,invading_fluid=water2,defending_fluid=air2,inlets=inlets,outlets=outlets)
#
##----------------------------------------------------------------------
#'''Performm a Diffusion Simulation on Partially Filled Network'''
##----------------------------------------------------------------------
##Apply desired/necessary pore scale physics methods
#air.regenerate()
#water.regenerate()
#OpenPNM.Physics.multi_phase.update_occupancy_OP(water,Pc=8000)
#OpenPNM.Physics.multi_phase.effective_occupancy(pn,air)
#OpenPNM.Physics.multi_phase.DiffusiveConductance(pn,air)
##Initialize algorithm object
#Fickian_alg = OpenPNM.Algorithms.FickianDiffusion()
##Create boundary condition arrays
#BCtypes = sp.zeros(pn.get_num_pores())
#BCvalues = sp.zeros(pn.get_num_pores())
##Specify Dirichlet-type and assign values
#OP_1.update()
#Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg',network=pn)
#Fickian_alg.set_pore_info(prop='Dirichlet',locations=pn.get_pore_indices(subdomain=['top','bottom']),is_indices=True)
#Dir_pores = sp.zeros_like(pn.get_pore_indices(subdomain='top'))
#Dir_pores[pn.get_pore_indices(subdomain='top')] = 0.8
##Dir_pores[pn.get_pore_indices(subdomain='bottom')] = 0.2
#Fickian_alg.set_pore_data(subdomain='Dirichlet',prop='BCval',data=Dir_pores,indices=pn.get_pore_indices(subdomain=['top']))
##Neumann
##BCtypes[pn.pore_properties['type']==1] = 1
##BCtypes[pn.pore_properties['type']==6] = 4
##BCvalues[pn.pore_properties['type']==1] = 8e-1
##BCvalues[pn.pore_properties['type']==6] = 2e-10
#Fickian_alg.set_boundary_conditions(types=BCtypes,values=BCvalues)
##Run simulation
#Fickian_alg.run(active_fluid=air)
#
#
##Export to VTK
#OpenPNM.Visualization.VTK().write(pn,fluid=water)

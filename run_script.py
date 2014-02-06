import OpenPNM

#======================================================================
'''Build Topological Network'''
#======================================================================
pn = OpenPNM.Network.Cubic(name='cubic_1').generate(divisions=[35,35,35],lattice_spacing=[0.0001])

#======================================================================
'''Build Geometry'''
#======================================================================
geom_recipe = {
'name': 'stick_and_ball',
'pore_seed': {'method': 'random'},
'throat_seed': {'method': 'neighbor_min'},
'pore_diameter': {'method': 'sphere',
                 'name': 'weibull_min',
                 'shape': 2.5,
                 'loc': 6e-6,
                 'scale': 2e-5},
'throat_diameter': {'method': 'cylinder',
                   'name': 'weibull_min',
                   'shape': 2.5,
                   'loc': 6e-6,
                   'scale': 2e-5},
'pore_volume': {'method': 'sphere'},
'throat_volume': {'method': 'cylinder'},
'throat_length': {'method': 'straight'},
}
geom = OpenPNM.Geometry.GenericGeometry().create(network=pn,**geom_recipe)
#geom = OpenPNM.Geometry.GenericGeometry(network=pn,name='stick_and_ball')
#geom.add_method(prop='pore_seed',model='random')
#geom.add_method(prop='throat_seed',model='neighbor_min')
#geom.add_method(prop='pore_diameter',model='sphere',name='weibull_min',shape=2.5,loc='6e-6',scale=2e-5)
#geom.add_method(prop='throat_diameter',model='cylinder',name='weibull_min',shape=2.5,loc='6e-6',scale=2e-5)
#geom.add_method(prop='pore_volume',model='sphere')
#geom.add_method(prop='throat_volume',model='cylinder')
#geom.add_method(prop='throat_length',model='straight')

#======================================================================
'''Build Fluids'''
#======================================================================
#Define the fluids properties
air_recipe = {
'name': 'air',
'Pc': 3.771e6,
'Tc': 132.65,
'MW': 0.0291,
'diffusivity': {'method': 'Fuller',
                'MA': 0.03199,
                'MB': 0.0291,
                'vA': 16.3,
                'vB': 19.7},
'viscosity': {'method': 'Reynolds',
              'uo': 0.001,
              'b': 0.1},
'molar_density': {'method': 'ideal_gas',
                  'R': 8.314},
}
air = OpenPNM.Fluids.GenericFluid(loggername='AIR',loglevel=10).create(network=pn,**air_recipe)
#air = OpenPNM.Fluids.GenericFluid(loggername='AIR',loglevel=10,network=pn,name='air')
#air.set_pore_data(prop='Pc',data=132.65)
#air.set_pore_data(prop='Tc',data=3.771e6)
#air.set_pore_data(prop='MW',data=0.0291)
#air.add_method(prop='diffusivity',model='Fuller',MA=0.03199,MB=0.0291,vA=16.3,vB=19.7)
#air.add_method(prop='viscosity',model='Reynolds',uo=0.001,b=0.1)
#air.add_method(prop='molar_density',model='ideal_gas',R=8.314)

water_recipe = {
'name': 'water',
'Pc': 2.206e6,
'Tc': 647,
'MW': 0.0181,
'diffusivity': {'method': 'constant',
                'value': 1e-12},
'viscosity': {'method': 'constant',
              'value': 0.001},
'molar_density': {'method': 'constant',
                  'value': 44445},
'surface_tension': {'method': 'constant',
                    'value': 0.072},
'contact_angle': {'method': 'constant',
                  'value': 110},
}
water = OpenPNM.Fluids.GenericFluid(loggername='WATER',loglevel=10).create(network=pn,**water_recipe)
#water = OpenPNM.Fluids.GenericFluid(loggername='AIR',loglevel=10,network=pn,name='water')
#water.set_pore_data(prop='Pc',data=132.65)
#water.set_pore_data(prop='Tc',data=3.771e6)
#water.set_pore_data(prop='MW',data=0.0291)
#water.add_method(prop='diffusivity',model='constant',value=1e-12)
#water.add_method(prop='viscosity',model='constant',value=0.001)
#water.add_method(prop='molar_density',model='constant',value=44445)
#water.add_method(prop='surface_tension',model='constant',value=0.072)
#water.add_method(prop='contact_angle',model='constant',value=110)


#======================================================================
'''Build Physics Objects'''
#======================================================================

phys_water = OpenPNM.Physics.GenericPhysics(network=pn,fluid=water,name='standard_water_physics')
phys_water.add_method(prop='capillary_pressure',model='purcell',r_torioid='1.e-5')
phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_water.add_method(prop='diffusive_conductance',model='bulk_diffusion')

phys_air = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air,name='standard_air_physics')
phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
phys_air.add_method(prop='diffusive_conductance',model='bulk_diffusion')

#======================================================================
'''Begin Simulations'''
#======================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#----------------------------------------------------------------------
#Initialize algorithm object
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(loglevel=10,loggername="OP",name='OP_1')
a = pn.get_pore_indices(subdomain='bottom')
#Run algorithm
#OP_1.run(network=pn,invading_fluid='water',defending_fluid='air',inlets=a,npts=20,AL=True)

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
#Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(name='Fickian_alg',network=pn)
#Fickian_alg.set_pore_info(prop='Dirichlet',data=pn.get_pore_indices(subdomain=['top','bottom']),indices=True)
#Dir_pores = sp.zeros(pn.get_num_pores())
#Dir_pores[pn.get_pore_indices(subdomain='top')] = 0.8
#Dir_pores[pn.get_pore_indices(subdomain='bottom')] = 0.2
#Fickian_alg.set_pore_data(subdomain='Dirichlet',prop='BCval',data=Dir_pores)
##Neumann
##BCtypes[pn.pore_properties['type']==1] = 1
##BCtypes[pn.pore_properties['type']==6] = 4
##BCvalues[pn.pore_properties['type']==1] = 8e-1
##BCvalues[pn.pore_properties['type']==6] = 2e-10
#Fickian_alg.set_boundary_conditions(types=BCtypes,values=BCvalues)
##Run simulation
#Fickian_alg.run(pn,active_fluid=air)
#
#
##Export to VTK
#OpenPNM.Visualization.VTK().write(pn,fluid=water)

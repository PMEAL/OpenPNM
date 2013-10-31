
import OpenPNM
import scipy as sp
from time import clock
start=clock()

#======================================================================
'''Generate Network Geometry'''
#======================================================================
#Define generation parameters
params_geo1= {
'domain_size': [],  #physical network size [meters]
'divisions': [10,10,10], #Number of pores in each direction
'lattice_spacing': [0.1],  #spacing between pores [meters]
'stats_pores': {'name': 'weibull_min', #Each statistical package takes different params, so send as dict
                'shape': 1.5,
                'loc': 6e-6,
                'scale': 2e-5},
'stats_throats': {'name': 'weibull_min',
                  'shape': 1.5,
                  'loc': 6e-6,
                  'scale': 2e-5},
                  'btype': [0,0,0]  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
}
#Generate Geometry
pn = OpenPNM.Geometry.Cubic(loglevel=40).generate(**params_geo1)

#======================================================================
'''Generate Fluids'''
#======================================================================
#Define the fluids properties
air_recipe = {
'name': 'air',
'Pc': 3.771e6, #Pa
'Tc': 132.65,  #K
'MW': 0.0291,  #kg/mol
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
'surface_tension': {'method': 'constant',
                    'value': 0},
'contact_angle': {'method': 'na'},
}
water_recipe = {
'name': 'water',
'Pc': 2.206e6, #Pa
'Tc': 647,     #K
'MW': 0.0181,  #kg/mol
'diffusivity': {'method': 'constant',
                'value': 1e-12},
'viscosity': {'method': 'constant',
              'value': 0.001},
'molar_density': {'method': 'constant',
                  'value': 44445},
'surface_tension': {'method': 'Eotvos',
                    'k': 2.25e-4},
'contact_angle': {'method': 'constant',
                  'value': 120},
}
#Create fluids
air = OpenPNM.Fluids.GenericFluid(loglevel=50).create(air_recipe)
water= OpenPNM.Fluids.GenericFluid(loglevel=50).create(water_recipe)
#set water and air as a fluid pair
water.set_pair(air)
#Set desired base conditions in the Fluids
air.pore_conditions['temperature'] = 353
air.pore_conditions['pressure'] = 101325
water.pore_conditions['temperature'] = 353
water.pore_conditions['pressure'] = 101325
#Update Fluids to the new conditions
water.regenerate()
air.regenerate()

#======================================================================
'''Begin Simulations'''
#======================================================================
'''Peform a Drainage Experiment (OrdinaryPercolation)'''
#----------------------------------------------------------------------
#Initialize algorithm object
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation()
#Apply desired/necessary pore scale physics methods
OpenPNM.Physics.CapillaryPressure.Washburn(pn,water)
a = pn.pore_properties['type']==2
#Run algorithm
OP_1.run(network=pn,invading_fluid=water,defending_fluid=air,inlets=a,npts=50,AL=True)

b = pn.pore_properties['type']==5
OP_1.evaluate_trapping(network=pn,invading_fluid=water,outlets=b)

#----------------------------------------------------------------------
'''Perform an Injection Experiment (InvasionPercolation)'''
#----------------------------------------------------------------------
#Create some new Fluids
water2 = OpenPNM.Fluids.GenericFluid(loglevel=50).create(water_recipe)
air2 = OpenPNM.Fluids.GenericFluid(loglevel=50).create(air_recipe)
water.pore_conditions['temperature'] = 353
water.pore_conditions['pressure'] = 101325
air.pore_conditions['temperature'] = 353
air.pore_conditions['pressure'] = 101325
#Initialize algorithm object
IP_1 = OpenPNM.Algorithms.InvasionPercolation()
#Apply desired/necessary pore scale physics methods
OpenPNM.Physics.CapillaryPressure.Washburn(pn,water2)
face = pn.pore_properties['type']==3
quarter = sp.rand(pn.get_num_pores(),)<.1
inlets = pn.pore_properties['numbering'][face&quarter]
outlets = pn.pore_properties['numbering'][pn.pore_properties['type']==4]
IP_1.run(pn,invading_fluid=water2,defending_fluid=air2,inlets=inlets,outlets=outlets)

#----------------------------------------------------------------------
'''Performm a Diffusion Simulation on Partially Filled Network'''
#----------------------------------------------------------------------
#Apply desired/necessary pore scale physics methods
air.regenerate()
water.regenerate()
OpenPNM.Physics.MultiPhase.update_occupancy_OP(water,Pc=8000)
OpenPNM.Physics.MultiPhase.effective_occupancy(pn,air)
OpenPNM.Physics.MassTransport.DiffusiveConductance(pn,air)
#Initialize algorithm object
Fickian_alg = OpenPNM.Algorithms.FickianDiffusion()
#Create boundary condition arrays
BCtypes = sp.zeros(pn.get_num_pores())
BCvalues = sp.zeros(pn.get_num_pores())
#Specify Dirichlet-type and assign values
BCtypes[pn.pore_properties['type']==1] = 1
BCtypes[pn.pore_properties['type']==6] = 1
BCvalues[pn.pore_properties['type']==1] = 8e-2
BCvalues[pn.pore_properties['type']==6] = 8e-1
#Neumann
#BCtypes[pn.pore_properties['type']==1] = 1
#BCtypes[pn.pore_properties['type']==6] = 4
#BCvalues[pn.pore_properties['type']==1] = 8e-1
#BCvalues[pn.pore_properties['type']==6] = 2e-10
Fickian_alg.set_boundary_conditions(types=BCtypes,values=BCvalues)
#Run simulation
Fickian_alg.run(pn,active_fluid=air)


#Export to VTK
OpenPNM.Visualization.VTK().write(pn,fluid=water)

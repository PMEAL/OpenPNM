import OpenPNM
import scipy as sp

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(name='cubic_1',loglevel=0).generate(divisions=[20, 20, 20], lattice_spacing=[0.0001],add_boundaries=True)

#==============================================================================
'''Build Geometry'''
#==============================================================================
geom = OpenPNM.Geometry.Stick_and_Ball(network=pn, name='geom', loglevel=0)
geom.regenerate()

#==============================================================================
'''Build Fluids'''
#==============================================================================
air = OpenPNM.Fluids.Air(network=pn, loglevel=0)
air.apply_conditions(temperature=350, pressure=200000)
air.regenerate()

water = OpenPNM.Fluids.Water(network=pn,loglevel=0)
water.add_method(prop='diffusivity',prop_name='DAB',model='constant',value=5e-12)
water.regenerate()

#==============================================================================
'''Build Physics Objects'''
#==============================================================================
phys_water = OpenPNM.Physics.GenericPhysics(network=pn, fluid=water,geometry=geom,name='phys_water', loglevel=0)
phys_water.add_method(prop='capillary_pressure', model='washburn')
phys_water.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_water.add_method(prop='diffusive_conductance', prop_name='gdAB', model='bulk_diffusion', diffusivity='DAB')
phys_water.regenerate()

phys_air = OpenPNM.Physics.GenericPhysics(network=pn, fluid=air,geometry=geom, name='phys_air', loglevel=0)
phys_air.add_method(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_air.add_method(prop='diffusive_conductance', model='bulk_diffusion')
phys_air.regenerate()

#==============================================================================
'''Begin Simulations'''
#==============================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#------------------------------------------------------------------------------
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(loglevel=0,loggername='OP',name='OP_1',network=pn)
a = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
OP_1.setup(invading_fluid='water',defending_fluid='air',inlets=a,npts=20)
OP_1.run()
#OP_1.plot_drainage_curve()

#------------------------------------------------------------------------------
'''Perform Fickian Diffusion'''
#------------------------------------------------------------------------------
Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loglevel=0, loggername='Fickian', name='Fickian_alg',network=pn)
# Assign Dirichlet boundary conditions to top and bottom surface pores
BC1_pores = pn.get_pore_indices(labels=['top','boundary'],mode='intersection')
Fickian_alg.set_pore_info(label='Dirichlet', locations=BC1_pores)
Fickian_alg.set_pore_data(prop='BCval', data=0.6, locations=BC1_pores)
BC2_pores = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
Fickian_alg.set_pore_info(label='Dirichlet', locations=BC2_pores)
Fickian_alg.set_pore_data(prop='BCval', data=0.2, locations=BC2_pores)

# Updating data based on the result of Percolation Algorithms
OP_1.update(Pc=3000)
# Run simulation
Fickian_alg.run(active_fluid=air)
Fickian_alg.update()

#------------------------------------------------------------------------------
#Export to VTK
import time, os

def timeit(fn, *args, **kwargs):
    start = time.time()
    output = fn(*args, **kwargs)
    print "{:.5f}s".format( time.time() - start )
    return output

timeit(OpenPNM.Visualization.VTK().write, net=pn, fluids=[air,water], filename='test.vtp') # 1.28s
timeit(OpenPNM.Visualization.Vtp.write, network=pn, fluids=[air, water], filename='test2.vtp') # 0.89s
print "These should be similar"
print os.stat('test.vtp').st_size
print os.stat('test2.vtp').st_size
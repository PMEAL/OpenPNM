import OpenPNM
import scipy as sp

#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(name='cubic_1',loglevel=20)
pn.generate(divisions=[20, 20, 20], lattice_spacing=[0.0001],add_boundaries=True)

#==============================================================================
'''Build Geometry'''
#==============================================================================
geom = OpenPNM.Geometry.Toray090(network=pn)
geom.set_locations(pores='all',throats='all')
geom.regenerate()

#==============================================================================
'''Build Fluids'''
#==============================================================================
air = OpenPNM.Fluids.Air(network=pn, loglevel=20)
air.apply_conditions(temperature=350, pressure=200000)
air.add_property(prop='electrical_conductivity',model='constant',value=5e-12)
air.regenerate()

water = OpenPNM.Fluids.Water(network=pn,loglevel=20)
water.add_property(prop='diffusivity',prop_name='DAB',model='constant',value=5e-12)
water.regenerate()

#==============================================================================
'''Build Physics Objects'''
#==============================================================================
phys_water = OpenPNM.Physics.GenericPhysics(network=pn, fluid=water,geometry=geom,name='phys_water')
phys_water.add_property(prop='capillary_pressure', model='washburn')
phys_water.add_property(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_water.add_property(prop='diffusive_conductance', prop_name='gdAB', model='bulk_diffusion', diffusivity='DAB')
phys_water.regenerate()

phys_air = OpenPNM.Physics.GenericPhysics(network=pn, fluid=air,geometry=geom, name='phys_air')
phys_air.add_property(prop='hydraulic_conductance', model='hagen_poiseuille')
phys_air.add_property(prop='diffusive_conductance', model='bulk_diffusion')
phys_air.add_property(prop='electronic_conductance', model='series_resistors')
phys_air.regenerate()

#==============================================================================
'''Begin Simulations'''
#==============================================================================
'''Perform a Drainage Experiment (OrdinaryPercolation)'''
#------------------------------------------------------------------------------
OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(loglevel=20,loggername='OP',name='OP_1',network=pn)
a = pn.pores(labels=['bottom','boundary'],mode='intersection')
OP_1.setup(invading_fluid=water,defending_fluid=air,inlets=a,npts=20)
OP_1.run()
#OP_1.plot_drainage_curve()

#------------------------------------------------------------------------------
'''Perform Fickian Diffusion'''
#------------------------------------------------------------------------------
Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loglevel=20, loggername='Fickian', name='Fickian_alg',network=pn)
# Assign Dirichlet boundary conditions to top and bottom surface pores
BC1_pores = pn.pores(labels=['top','front'],mode='intersection')
Fickian_alg.set_info(label='Dirichlet', pores=BC1_pores)
Fickian_alg.set_data(prop='BCval', data=0.6, pores=BC1_pores)
BC2_pores = pn.pores(labels=['top','back'],mode='intersection')
Fickian_alg.set_info(label='Dirichlet', pores=BC2_pores)
Fickian_alg.set_data(prop='BCval', data=0.2, pores=BC2_pores)

# Updating data based on the result of Percolation Algorithms
OP_1.update(Pc=11000)
# Run simulation
Fickian_alg.run(active_fluid=air)
Fickian_alg.update()

#------------------------------------------------------------------------------
'''Export to VTK'''
#------------------------------------------------------------------------------
#OpenPNM.Visualization.VTK().write(net=pn, fluids=[air,water], filename='output.vtp')
OpenPNM.Visualization.Vtp.write(filename='test.vtp',fluids=[air,water],network=pn)










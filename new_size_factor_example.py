import openpnm as op
import numpy as np

pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4)
geo = op.geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
geo.settings['interpolation_mode'] = 'min'

geo.add_model(propname='pore.seed',
              element='pore',
              model=op.models.misc.random,
              num_range=[0.1, 0.9])

geo.add_model(propname='pore.diameter',
              model=op.models.geometry.pore_size.weibull,
              shape=1.5, scale=1e-5, loc=5e-5)

geo.add_model(propname='pore.area',
              model=op.models.geometry.pore_area.sphere)

geo.add_model(propname='throat.inherited_diameter',
              model=op.models.geometry.throat_size.from_neighbor_pores,
              mode='min')

geo.add_model(propname='throat.diameter',
              model=op.models.misc.scaled,
              prop='throat.inherited_diameter',
              factor=0.5)


geo.add_model(propname='throat.area',
              model=op.models.geometry.throat_area.cylinder)

geo.add_model(propname='throat.endpoints',
              model=op.models.geometry.throat_endpoints.spherical_pores)

geo.add_model(propname='throat.conduit_length',
              model=op.models.geometry.throat_length.conduit_lengths)

mod_g = op.models.geometry.conduit_hydraulic_coefficient.spheres_and_cylinders
geo.add_model(propname='throat.flow_coeff',
              model=mod_g,
              conduit_lengths='throat.conduit_length')

air = op.phases.Air(network=pn, name='air')

# new model with conduit length
phys = op.physics.GenericPhysics(network=pn, phase=air, geometry=geo)
mod = op.models.physics.hydraulic_conductance.generic_hydraulic
phys.add_model(propname='throat.hydraulic_conductance_new',
               model=mod, shape_coeff='throat.flow_coeff')
hc_new_w_cl=phys['throat.hydraulic_conductance_new']

# new model without conduit length

geo.add_model(propname='throat.flow_coeff_wo_cl',
              model=mod_g)
phys.add_model(propname='throat.hydraulic_conductance_new_wo_cl',
               model=mod, shape_coeff='throat.flow_coeff_wo_cl')
hc_new_wo_cl=phys['throat.hydraulic_conductance_new_wo_cl']

# old model using shape factors
mod = op.models.geometry.hydraulic_shape_factors.ball_and_stick
geo.add_model(propname='throat.flow_shape_factor',
              model=mod)
mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys.add_model(propname='throat.hydraulic_conductance_old',
               model=mod, conduit_shape_factors='throat.flow_shape_factor',
               conduit_lengths='throat.conduit_length')
hc_old=phys['throat.hydraulic_conductance_old']

print('avg old model '+ str(np.mean(hc_old))+ 'new model with CL'+ str(np.mean(hc_new_w_cl))+ 'without cl'+ str(np.mean(hc_new_wo_cl)) )




# sf = op.algorithms.StokesFlow(network=pn, phase=air)
# sf.setup(conductance='throat.hydraulic_conductance_new')
# sf.set_value_BC(pores=pn.pores('left'), values=101323)
# sf.set_value_BC(pores=pn.pores('right'), values=0)
# sf.run()
# Q1 = sf.rate(pores=pn.pores('left'))
# print(Q1)

# sf = op.algorithms.StokesFlow(network=pn, phase=air)
# sf.setup(conductance='throat.hydraulic_conductance_new')
# sf.set_value_BC(pores=pn.pores('left'), values=101323)
# sf.set_value_BC(pores=pn.pores('right'), values=0)
# sf.run()
# Q2 = sf.rate(pores=pn.pores('left'))
# print(Q2)

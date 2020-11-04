# import openpnm as op
# import numpy as np

# pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4)
# geo = op.geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
# geo.settings['interpolation_mode'] = 'min'

# geo.add_model(propname='pore.seed',
#               element='pore',
#               model=op.models.misc.random,
#               num_range=[0.1, 0.9])

# geo.add_model(propname='pore.diameter',
#               model=op.models.geometry.pore_size.weibull,
#               shape=1.5, scale=1e-5, loc=5e-5)

# geo.add_model(propname='throat.inherited_diameter',
#               model=op.models.geometry.throat_size.from_neighbor_pores,
#               mode='min')

# geo.add_model(propname='throat.diameter',
#               model=op.models.misc.scaled,
#               prop='throat.inherited_diameter',
#               factor=0.5)

# geo.add_model(propname='throat.endpoints',
#               model=op.models.geometry.throat_endpoints.spherical_pores)

# geo.add_model(propname='throat.conduit_length',
#               model=op.models.geometry.throat_length.conduit_lengths)

# mod = op.models.geometry.conduit_diffusive_coefficient.spheres_and_cylinders
# geo.add_model(propname='throat.diff_coeff_new_wo_cl',
#               model=mod,
#               conduit_lengths=None)
# geo.add_model(propname='throat.diff_coeff_new_w_cl',
#               model=mod,
#               conduit_lengths='throat.conduit_length')
# print('the difference is'+ str(np.mean(geo['throat.diff_coeff_new_wo_cl']-geo['throat.diff_coeff_new_w_cl'])))
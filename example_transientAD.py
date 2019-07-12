import openpnm as op
import numpy as np
ws = op.Workspace()
proj = ws.new_project()
# ws.settings['loglevel'] = 20


# network, geometry, phase
np.random.seed(0)

net = op.network.Cubic(shape=[8, 8, 1], spacing=9e-4, project=proj)
prs = (net['pore.back'] * net['pore.right'] + net['pore.back'] *
       net['pore.left'] + net['pore.front'] * net['pore.right'] +
       net['pore.front'] * net['pore.left'])
thrts = net['throat.surface']
op.topotools.trim(network=net, pores=net.Ps[prs], throats=net.Ts[thrts])


geo = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)
pore_d = op.models.misc.constant
throat_d = op.models.misc.constant
geo.add_model(propname='pore.diameter', model=pore_d, value=1.5e-4)
geo.add_model(propname='throat.diameter', model=throat_d, value=1e-4)
geo.regenerate_models()

water = op.phases.Water(network=net)


# physics
phys = op.physics.GenericPhysics(network=net, phase=water, geometry=geo)

flow = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys.add_model(propname='throat.hydraulic_conductance',
               pore_viscosity='pore.viscosity',
               throat_viscosity='throat.viscosity',
               model=flow, regen_mode='normal')

# algorithms
sf = op.algorithms.StokesFlow(network=net, phase=water)
sf.set_value_BC(pores=net.pores('back'), values=0.01)
sf.set_value_BC(pores=net.pores('front'), values=0.00)
sf.settings['rxn_tolerance'] = 1e-12
sf.run()
water.update(sf.results())

dif = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance', model=dif,
               regen_mode='normal')

ad_dif = op.models.physics.ad_dif_conductance.ad_dif
phys.add_model(propname='throat.ad_dif_conductance',
               model=ad_dif, regen_mode='normal')

ad = op.algorithms.TransientAdvectionDiffusion(network=net, phase=water)
ad.set_value_BC(pores=net.pores('back'), values=100)
ad.set_value_BC(pores=net.pores('front'), values=90)
ad.settings['rxn_tolerance'] = 1e-12

ad.settings['i_max_iter'] = 10
ad.settings['i_tolerance'] = 1e-04
ad.settings['charge_conservation'] = 'laplace'
ad.settings['t_output'] = 500
ad.settings['t_step'] = 100
ad.settings['t_final'] = 2000
# Electroneutrality condition does not work with new Mixtures
# pnp.settings['charge_conservation'] = 'electroneutrality'
ad.run()

water.update(ad.results())

# output results to a vtk file
# proj.export_data(phases=[water], filename='OUT', filetype='xdmf')

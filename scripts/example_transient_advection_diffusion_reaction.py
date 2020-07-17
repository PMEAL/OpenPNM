import openpnm as op
import numpy as np
ws = op.Workspace()
proj = ws.new_project()
# ws.settings['loglevel'] = 20


"""
    Details about the continum and numerical model equations can be found on:
    Sadeghi, M. A., Agnaou, M., Barralet, J., & Gostick, J. (2020).
    Dispersion modeling in pore networks: A comparison of common pore-scale
    models and alternative approaches.
    Journal of Contaminant Hydrology, 228, 103578.
"""


# network
np.random.seed(0)
net = op.network.Cubic(shape=[33, 33, 1], spacing=9e-4, project=proj)
# remove the pores at the corners of the network
prs = (net['pore.back'] * net['pore.right'] + net['pore.back']
       * net['pore.left'] + net['pore.front'] * net['pore.right']
       + net['pore.front'] * net['pore.left'])
thrts = net['throat.surface']
op.topotools.trim(network=net, pores=net.Ps[prs], throats=net.Ts[thrts])


# geometry
geo = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)
# define constant pore and throat diameters
pore_d = op.models.misc.constant
throat_d = op.models.misc.constant
geo.add_model(propname='pore.diameter', model=pore_d, value=1.5e-4)
geo.add_model(propname='throat.diameter', model=throat_d, value=1e-4)
geo.regenerate_models()


# phase
water = op.phases.Water(network=net)


# physics
phys = op.physics.GenericPhysics(network=net, phase=water, geometry=geo)
# Stokes flow
flow = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys.add_model(propname='throat.hydraulic_conductance',
               pore_viscosity='pore.viscosity',
               throat_viscosity='throat.viscosity',
               model=flow, regen_mode='normal')
# diffusion
dif = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance', model=dif,
               regen_mode='normal')
# advection diffusion
ad_dif = op.models.physics.ad_dif_conductance.ad_dif
phys.add_model(propname='throat.ad_dif_conductance',
               model=ad_dif, regen_mode='normal')
# source term
linear = op.models.physics.generic_source_term.linear
phys['pore.A1'] = -1e-15
phys['pore.A2'] = 0.0
phys.add_model(propname='pore.rxn', model=linear, X='pore.concentration',
               A1='pore.A1', A2='pore.A2')
rxn_pores = np.array([200])
net.set_label('rxn', pores=rxn_pores)


# algorithms
# Stokes flow algorithm
sf = op.algorithms.StokesFlow(network=net, phase=water)
sf.set_value_BC(pores=net.pores('back'), values=0.01)
sf.set_value_BC(pores=net.pores('front'), values=0.00)
sf.run()
water.update(sf.results())
# transient advection diffusion
ad = op.algorithms.TransientAdvectionDiffusion(network=net, phase=water)
ad.set_source(propname='pore.rxn', pores=rxn_pores)
ad.set_value_BC(pores=net.pores('back'), values=100)
ad.set_value_BC(pores=net.pores('front'), values=90)
# change 'steady' to 'implicit' or 'cranknicolson' for transient simulations
settings = {'t_scheme': 'steady', 'solver_tol': 1e-12, 't_output': 1000,
            't_step': 500, 't_final': 100000}
ad.settings.update(settings)
ad.set_IC(90)
ad.run()
water.update(ad.results())


# output results for visualization on Paraview
# proj.export_data(phases=[water], filename='OUT', filetype='xdmf')

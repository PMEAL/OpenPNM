import openpnm as op
from openpnm.phases import mixtures
import numpy as np
ws = op.Workspace()
proj = ws.new_project()
# ws.settings['loglevel'] = 20


"""
    Details about the continum and numerical model equations can be found on:
    Agnaou, M., Sadeghi, M. A., Tranter, T. G., & Gostick, J. (2020).
    Modeling transport of charged species in pore networks: solution of the
    Nernst-Planck equations coupled with fluid flow and charge conservation
    equations.
    Computers & Geosciences, 104505.
"""


# network, geometry, phase
np.random.seed(0)

net = op.network.Cubic(shape=[23, 15, 1], spacing=1e-6, project=proj)
prs = (net['pore.back'] * net['pore.right'] + net['pore.back']
       * net['pore.left'] + net['pore.front'] * net['pore.right']
       + net['pore.front'] * net['pore.left'])
prs = net.Ps[prs]

thrts = net['throat.surface']
thrts = net.Ts[thrts]

op.topotools.trim(network=net, pores=prs, throats=thrts)

np.random.seed(0)
op.topotools.reduce_coordination(net, 3)

np.random.seed(0)
geo = op.geometry._StickAndBall(network=net, pores=net.Ps, throats=net.Ts)


sw = mixtures.SalineWater(network=net)
# Retrieve handles to each species for use below
Na = sw.components['Na_' + sw.name]
Cl = sw.components['Cl_' + sw.name]
H2O = sw.components['H2O_' + sw.name]

# physics
phys = op.physics.GenericPhysics(network=net, phase=sw, geometry=geo)

flow = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys.add_model(propname='throat.hydraulic_conductance',
               pore_viscosity='pore.viscosity',
               throat_viscosity='throat.viscosity',
               model=flow, regen_mode='normal')

current = op.models.physics.ionic_conductance.poisson
phys.add_model(propname='throat.ionic_conductance',
               model=current, regen_mode='normal')

eA_dif = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance.' + Na.name,
               pore_diffusivity='pore.diffusivity.' + Na.name,
               throat_diffusivity='throat.diffusivity.' + Na.name,
               model=eA_dif, regen_mode='normal')

eB_dif = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance.' + Cl.name,
               pore_diffusivity='pore.diffusivity.' + Cl.name,
               throat_diffusivity='throat.diffusivity.' + Cl.name,
               model=eB_dif, regen_mode='normal')

scheme = 'powerlaw'
ad_dif_mig_Na = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.' + Na.name,
               pore_pressure='pore.pressure', model=ad_dif_mig_Na,
               ion=Na.name, s_scheme=scheme)

ad_dif_Na = op.models.physics.ad_dif_conductance.ad_dif
phys.add_model(propname='throat.ad_dif_conductance.' + Na.name,
               throat_diffusive_conductance='throat.diffusive_conductance.' +
               Na.name, pore_pressure='pore.pressure', model=ad_dif_Na,
               s_scheme=scheme)

ad_dif_mig_Cl = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.' + Cl.name,
               pore_pressure='pore.pressure', model=ad_dif_mig_Cl,
               ion=Cl.name, s_scheme=scheme)

# settings for algorithms
setts1 = {'solver_family': 'scipy', 'solver_max_iter': 5, 'solver_tol': 1e-08,
          'solver_rtol': 1e-08, 'nlin_max_iter': 10, 'cache_A': False,
          'cache_b': False}
setts2 = {'g_tol': 1e-4, 'g_max_iter': 100}

# algorithms
sf = op.algorithms.StokesFlow(network=net, phase=sw, settings=setts1)
sf.set_value_BC(pores=net.pores('back'), values=6)
sf.set_value_BC(pores=net.pores('front'), values=1)
sf.run()
sw.update(sf.results())

p = op.algorithms.IonicConduction(network=net, phase=sw, settings=setts1)
p.set_value_BC(pores=net.pores('left'), values=0.01)
p.set_value_BC(pores=net.pores('right'), values=0.01)
p.settings['charge_conservation'] = 'poisson'
p.run()
sw.update(p.results())

eA = op.algorithms.NernstPlanck(network=net, phase=sw, ion=Na.name,
                                settings=setts1)
eA.set_value_BC(pores=net.pores('back'), values=20)
eA.set_value_BC(pores=net.pores('front'), values=10)
eA.run()
sw['pore.val1'] = eA['pore.concentration.Na_mix_01']

eB = op.algorithms.NernstPlanck(network=net, phase=sw, ion=Na.name,
                                settings=setts1)
eB.set_value_BC(pores=net.pores('back'), values=1)
eB.set_value_BC(pores=net.pores('left'), values=0.1)
eB.set_value_BC(pores=net.pores('right'), values=0.1)
eB.set_outflow_BC(pores=net.pores('front'))
eB.run()
sw['pore.val2'] = eB['pore.concentration.Na_mix_01']


eC = op.algorithms.AdvectionDiffusion(network=net, phase=sw, settings=setts1)
setts3 = setts1.copy()
setts3.update(
        {'conductance': 'throat.ad_dif_conductance.' + Na.name,
         'diffusive_conductance': 'throat.diffusive_conductance.' + Na.name,
         'hydraulic_conductance': 'throat.hydraulic_conductance'})
eC.settings.update(setts3)
eC.set_value_BC(pores=net.pores('back'), values=1)
eC.set_value_BC(pores=net.pores('left'), values=0.1)
eC.set_value_BC(pores=net.pores('right'), values=0.1)
eC.set_outflow_BC(pores=net.pores('front'))
eC.run()
sw['pore.val3'] = eC['pore.concentration']

# output data to Paraview
# sw['pore.val2'] should be equal to sw['pore.val3']
# proj.export_data(phases=[sw], filename='OUT', filetype='vtp')

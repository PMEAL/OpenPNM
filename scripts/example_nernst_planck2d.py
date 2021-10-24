r"""
Details about the continum and numerical model equations can be found on:
Agnaou, M., Sadeghi, M. A., Tranter, T. G., & Gostick, J. (2020).

Modeling transport of charged species in pore networks: solution of the
Nernst-Planck equations coupled with fluid flow and charge conservation
equations. Computers & Geosciences, 104505.

"""
import openpnm as op
from openpnm.phases import mixtures
import numpy as np


ws = op.Workspace()
proj = ws.new_project()
export = False

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
geo = op.geometry.StickAndBall2D(network=net, pores=net.Ps, throats=net.Ts)


sw = mixtures.SalineWater(network=net)
# Retrieve handles to each species for use below
Na = sw.components['Na_' + sw.name]
Cl = sw.components['Cl_' + sw.name]
H2O = sw.components['H2O_' + sw.name]

# physics
phys = op.physics.GenericPhysics(network=net, phase=sw, geometry=geo)

flow = op.models.physics.hydraulic_conductance.hagen_poiseuille_2d
phys.add_model(propname='throat.hydraulic_conductance',
               pore_viscosity='pore.viscosity',
               throat_viscosity='throat.viscosity',
               model=flow, regen_mode='normal')

current = op.models.physics.ionic_conductance.electroneutrality
phys.add_model(propname='throat.ionic_conductance', ions=[Na.name, Cl.name],
               model=current, regen_mode='normal')

eA_dif = op.models.physics.diffusive_conductance.ordinary_diffusion_2d
phys.add_model(propname='throat.diffusive_conductance.' + Na.name,
               pore_diffusivity='pore.diffusivity.' + Na.name,
               throat_diffusivity='throat.diffusivity.' + Na.name,
               model=eA_dif, regen_mode='normal')

eB_dif = op.models.physics.diffusive_conductance.ordinary_diffusion_2d
phys.add_model(propname='throat.diffusive_conductance.' + Cl.name,
               pore_diffusivity='pore.diffusivity.' + Cl.name,
               throat_diffusivity='throat.diffusivity.' + Cl.name,
               model=eB_dif, regen_mode='normal')

scheme = 'powerlaw'
ad_dif_mig_Na = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.' + Na.name,
               pore_pressure='pore.pressure', model=ad_dif_mig_Na,
               ion=Na.name, s_scheme=scheme)

ad_dif_mig_Cl = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.' + Cl.name,
               pore_pressure='pore.pressure', model=ad_dif_mig_Cl,
               ion=Cl.name, s_scheme=scheme)

# settings for algorithms
setts1 = {'solver_max_iter': 5, 'solver_tol': 1e-08, 'solver_rtol': 1e-08,
          'nlin_max_iter': 10, 'cache_A': False, 'cache_b': False}
setts2 = {'g_tol': 1e-4, 'g_max_iter': 100}

# algorithms
sf = op.algorithms.StokesFlow(network=net, phase=sw, settings=setts1)
sf.set_value_BC(pores=net.pores('back'), values=11)
sf.set_value_BC(pores=net.pores('front'), values=10)
sf.run()
sw.update(sf.results())

p = op.algorithms.IonicConduction(network=net, phase=sw, settings=setts1)
p.set_value_BC(pores=net.pores('left'), values=0.02)
p.set_value_BC(pores=net.pores('right'), values=0.01)
p.settings['charge_conservation'] = 'electroneutrality_2D'

eA = op.algorithms.NernstPlanck(network=net, phase=sw, ion=Na.name,
                                settings=setts1)
eA.set_value_BC(pores=net.pores('back'), values=20)
eA.set_value_BC(pores=net.pores('front'), values=10)

eB = op.algorithms.NernstPlanck(network=net, phase=sw, ion=Cl.name,
                                settings=setts1)
eB.set_value_BC(pores=net.pores('back'), values=20)
eB.set_value_BC(pores=net.pores('front'), values=10)

pnp = op.algorithms.NernstPlanckMultiphysicsSolver(network=net, phase=sw,
                                                   settings=setts2)
pnp.settings['potential_field'] = p.name
pnp.settings['ions'] = [eA.name, eB.name]
pnp.run()

sw.update(sf.results())
sw.update(p.results())
sw.update(eA.results())
sw.update(eB.results())

# output data to Paraview
if export:
    proj.export_data(phases=[sw], filename='out', filetype='xdmf')

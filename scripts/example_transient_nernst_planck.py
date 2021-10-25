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

net = op.network.Cubic(shape=[8, 8, 1], spacing=9e-4, project=proj)
prs = (net['pore.back'] * net['pore.right'] + net['pore.back']
       * net['pore.left'] + net['pore.front'] * net['pore.right']
       + net['pore.front'] * net['pore.left'])
thrts = net['throat.surface']
op.topotools.trim(network=net, pores=net.Ps[prs], throats=net.Ts[thrts])


geo = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)
pore_d = op.models.misc.constant
throat_d = op.models.misc.constant
geo.add_model(propname='pore.diameter', model=pore_d, value=1.5e-4)
geo.add_model(propname='throat.diameter', model=throat_d, value=1e-4)
geo.regenerate_models()

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

current = op.models.physics.ionic_conductance.electroneutrality
phys.add_model(propname='throat.ionic_conductance', ions=[Na.name, Cl.name],
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

s_scheme = 'powerlaw'
ad_dif_mig_Na = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.' + Na.name,
               pore_pressure='pore.pressure', model=ad_dif_mig_Na,
               ion=Na.name, s_scheme=s_scheme)

ad_dif_mig_Cl = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.' + Cl.name,
               pore_pressure='pore.pressure', model=ad_dif_mig_Cl,
               ion=Cl.name, s_scheme=s_scheme)

# settings for algorithms
setts1 = {'solver_max_iter': 5, 'solver_tol': 1e-08, 'solver_rtol': 1e-08,
          'nlin_max_iter': 10, 'cache_A': False, 'cache_b': False}
setts2 = {'g_tol': 1e-4, 'g_max_iter': 4, 't_output': 5000, 't_step': 500,
          't_final': 20000, 't_scheme': 'implicit'}

# algorithms
sf = op.algorithms.StokesFlow(network=net, phase=sw, settings=setts1)
sf.set_value_BC(pores=net.pores('back'), values=0.01)
sf.set_value_BC(pores=net.pores('front'), values=0.00)
sf.run()
sw.update(sf.results())

p = op.algorithms.TransientIonicConduction(network=net, phase=sw,
                                           settings=setts1)
p.set_value_BC(pores=net.pores('left'), values=0.1)
p.set_value_BC(pores=net.pores('right'), values=0.00)
p.settings['charge_conservation'] = 'electroneutrality'

eA = op.algorithms.TransientNernstPlanck(network=net, phase=sw, ion=Na.name,
                                         settings=setts1)
eA.set_value_BC(pores=net.pores('back'), values=100)
eA.set_value_BC(pores=net.pores('front'), values=90)

eB = op.algorithms.TransientNernstPlanck(network=net, phase=sw, ion=Cl.name,
                                         settings=setts1)
eB.set_value_BC(pores=net.pores('back'), values=100)
eB.set_value_BC(pores=net.pores('front'), values=90)

it = op.algorithms.TransientNernstPlanckMultiphysicsSolver(network=net,
                                                           phase=sw,
                                                           settings=setts2)
it.setup(potential_field=p.name, ions=[eA.name, eB.name])
it.run()

sw.update(sf.results())
sw.update(p.results())
sw.update(eA.results())
sw.update(eB.results())

# output results to a vtk file for visualization on Paraview
if export:
    proj.export_data(phases=[sw], filename='out', filetype='xdmf')

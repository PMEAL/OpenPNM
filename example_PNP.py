import openpnm as op
import numpy as np
ws = op.Workspace()
proj = ws.new_project()


# network, geometry, phase
np.random.seed(0)

net = op.network.Cubic(shape=[8, 8, 1], spacing=9e-4)
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

sw = op.phases.Mixtures.SalineWater(network=net)

# physics
phys = op.physics.GenericPhysics(network=net, phase=sw, geometry=geo)

flow = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys.add_model(propname='throat.hydraulic_conductance',
               pore_viscosity='pore.viscosity',
               throat_viscosity='throat.viscosity',
               model=flow, regen_mode='normal')

current = op.models.physics.ionic_conductance.laplace
phys.add_model(propname='throat.ionic_conductance',
               model=current, regen_mode='normal')

eA_dif = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance.Na',
               pore_diffusivity='pore.diffusivity.Na',
               throat_diffusivity='throat.diffusivity.Na',
               model=eA_dif, regen_mode='normal')

eB_dif = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance.Cl',
               pore_diffusivity='pore.diffusivity.Cl',
               throat_diffusivity='throat.diffusivity.Cl',
               model=eB_dif, regen_mode='normal')

# algorithms
sf = op.algorithms.StokesFlow(network=net, phase=sw)
sf.set_value_BC(pores=net.pores('back'), values=0.01)
sf.set_value_BC(pores=net.pores('front'), values=0.00)
sf.settings['rxn_tolerance'] = 1e-12
sf.run()
sw.update(sf.results())

p = op.algorithms.OhmicConduction(network=net, phase=sw)
p.settings['conductance'] = 'throat.ionic_conductance'
p.settings['quantity'] = 'pore.potential'
p.set_value_BC(pores=net.pores('left'), values=0.01)
p.set_value_BC(pores=net.pores('right'), values=0.00)
p.settings['rxn_tolerance'] = 1e-12


eA = op.algorithms.NernstPlanck(network=net, phase=sw, ion='Na')
eA.set_value_BC(pores=net.pores('back'), values=100)
eA.set_value_BC(pores=net.pores('front'), values=90)
eA.settings['rxn_tolerance'] = 1e-12

eB = op.algorithms.NernstPlanck(network=net, phase=sw, ion='Cl')
eB.set_value_BC(pores=net.pores('back'), values=100)
eB.set_value_BC(pores=net.pores('front'), values=90)
eB.settings['rxn_tolerance'] = 1e-12

ad_dif_mig_Na = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.Na',
               model=ad_dif_mig_Na, ion='Na',
               s_scheme='exponential')

ad_dif_mig_Cl = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.Cl',
               pore_pressure=sf.settings['quantity'],
               model=ad_dif_mig_Cl, ion='Cl',
               s_scheme='exponential')

pnp = op.algorithms.ChargeConservationNernstPlanck(network=net, phase=sw)
pnp.setup(potential_field=p, ions=[eA, eB])
pnp.settings['max_iter'] = 10
pnp.settings['tolerance'] = 1e-04
pnp.settings['charge_conservation'] = 'laplace'
pnp.run()

sw.update(sf.results())
sw.update(p.results())
sw.update(eA.results())
sw.update(eB.results())

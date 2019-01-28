import openpnm as op
import numpy as np
ws = op.Workspace()
proj = ws.new_project()
# ws.settings['loglevel'] = 20


# network, geometry, phase
np.random.seed(0)
net = op.network.Cubic(shape=[13, 13, 1], spacing=1e-9)
geo = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)
sw = op.phases.Mixtures.SalineWater(network=net)

# physics
phys = op.physics.GenericPhysics(network=net, phase=sw, geometry=geo)

flow = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys.add_model(propname='throat.hydraulic_conductance',
               pore_viscosity='pore.viscosity',
               throat_viscosity='throat.viscosity',
               model=flow, regen_mode='normal')

current = op.models.physics.ionic_conductance.ordinary
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
sf.set_value_BC(pores=net.pores('back'), values=20)
sf.set_value_BC(pores=net.pores('front'), values=1.00)
sf.settings['rxn_tolerance'] = 1e-12
sw.update(sf.results())

p = op.algorithms.OhmicConduction(network=net, phase=sw)
p.settings['conductance'] = 'throat.ionic_conductance'
p.settings['quantity'] = 'pore.potential'
p.set_value_BC(pores=net.pores('left'), values=0.01)
p.set_value_BC(pores=net.pores('right'), values=0.0)
p.settings['rxn_tolerance'] = 1e-12
sw.update(p.results())

eA = op.algorithms.NernstPlanck(network=net, phase=sw, electrolyte='Na')
eA.set_value_BC(pores=net.pores('back'), values=300)
eA.set_value_BC(pores=net.pores('front'), values=100)
eA.settings['rxn_tolerance'] = 1e-12

eB = op.algorithms.NernstPlanck(network=net, phase=sw, electrolyte='Cl')
eB.set_value_BC(pores=net.pores('back'), values=300)
eB.set_value_BC(pores=net.pores('front'), values=100)
eB.settings['rxn_tolerance'] = 1e-12

ad_dif_mig_Na = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.Na',
               model=ad_dif_mig_Na, electrolyte='Na')

ad_dif_mig_Cl = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.Cl',
               model=ad_dif_mig_Cl, electrolyte='Cl')

pnp = op.algorithms.PoissonNernstPlanck(network=net, phase=sw)
pnp.setup(potential_field=p, electrolytes=[eA, eB])
pnp.settings['max_iter'] = 5
pnp.settings['tolerance'] = 1e-08
# pnp.settings['charge_conservation'] = 'poisson'
# pnp.run()

# sw.update(sf.results())
# sw.update(p.results())
# sw.update(eA.results())
# sw.update(eB.results())
# op.io.VTK.save(network=net, phases=[sw], filename='OUTPUT_PNP')

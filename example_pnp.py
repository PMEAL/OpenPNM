import openpnm as op
import numpy as np
ws = op.Workspace()
proj = ws.new_project()
ws.settings['loglevel'] = 20


# network, geometry, phase
np.random.seed(0)
net = op.network.Cubic(shape=[19, 13, 1], spacing=1e-7)
geo = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)
sw = op.phases.Mixtures.SalineWater(network=net)

# physics
phys = op.physics.GenericPhysics(network=net, phase=sw, geometry=geo)

flow = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys.add_model(propname='throat.hydraulic_conductance.solvent',
               pore_viscosity='pore.viscosity.solvent',
               throat_viscosity='throat.viscosity.solvent',
               model=flow, regen_mode='normal')
phys.regenerate_models()

current = op.models.physics.ionic_conductance.ordinary
phys.add_model(propname='throat.electrical_conductance.solvent',
               model=current, regen_mode='normal')
phys.regenerate_models()

eA_dif = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance.Na',
               pore_diffusivity='pore.diffusivity.Na',
               throat_diffusivity='throat.diffusivity.Na',
               model=eA_dif, regen_mode='normal')
phys.regenerate_models()

eB_dif = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance.Cl',
               pore_diffusivity='pore.diffusivity.Cl',
               throat_diffusivity='throat.diffusivity.Cl',
               model=eB_dif, regen_mode='normal')
phys.regenerate_models()

# algorithms
sf = op.algorithms.ReactiveTransport(network=net, phase=sw)
sf.set_value_BC(pores=net.pores('back'), values=0.2)
sf.set_value_BC(pores=net.pores('front'), values=0.1)
sf.settings['solver_family'] = 'pyamg'
sf.settings['rxn_tolerance'] = 1e-12

p = op.algorithms.ReactiveTransport(network=net, phase=sw)
p.set_value_BC(pores=net.pores('left'), values=0.1)
p.set_value_BC(pores=net.pores('right'), values=1.1)
p.settings['solver_family'] = 'pyamg'
p.settings['rxn_tolerance'] = 1e-12

eA = op.algorithms.ReactiveTransport(network=net, phase=sw, name='Na')
eA.set_value_BC(pores=net.pores('back'), values=5e-3)
eA.settings['rxn_tolerance'] = 1e-12

eB = op.algorithms.ReactiveTransport(network=net, phase=sw, name='Cl')
eB.set_value_BC(pores=net.pores('back'), values=5e-3)
eB.settings['rxn_tolerance'] = 1e-12

pnp = op.algorithms.PoissonNernstPlanck(network=net)
pnp.setup(phase=sw, pressure_field=sf,
          potential_field=p, electrolytes=[eA, eB])
pnp.settings['s_scheme'] = 'exponential'
pnp.settings['max_iter'] = 20
pnp.settings['tolerance'] = 1e-08
pnp.run()

sw.update(sf.results())
sw.update(p.results())
sw.update(eA.results())
sw.update(eB.results())
op.io.VTK.save(network=net, phases=[sw], filename='OUTPUT_PNP')

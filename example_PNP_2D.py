import openpnm as op
from openpnm.phases import mixtures
import numpy as np
ws = op.Workspace()
proj = ws.new_project()
# ws.settings['loglevel'] = 20


scheme = 'powerlaw'

# network, geometry, phase
np.random.seed(0)

net = op.network.Cubic(shape=[23, 15, 1], spacing=1e-6, project=proj)
prs = (net['pore.back'] * net['pore.right'] + net['pore.back'] *
       net['pore.left'] + net['pore.front'] * net['pore.right'] +
       net['pore.front'] * net['pore.left'])
prs = net.Ps[prs]

b_prs_1 = np.append(net.pores('back'), net.pores('front'))
b_prs_2 = np.append(net.pores('left'), net.pores('right'))
b_prs = np.append(b_prs_1, b_prs_2)
b_thrts = net.find_neighbor_throats(b_prs)

thrts_1 = net['throat.surface']
thrts_1 = net.Ts[thrts_1]
np.random.seed(0)
thrts_i = net.Ts[~net['throat.surface']]
thrts_sample = [i for i in thrts_i if i not in b_thrts]
L = int(0.05*len(thrts_sample))
thrts_2 = np.random.choice(thrts_sample, size=(L,), replace=False)
thrts_2 = np.array([])
thrts = np.append(thrts_1, thrts_2)

op.topotools.trim(network=net, pores=prs, throats=thrts)

op.topotools.reduce_coordination(net, 3)

geo = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)


sw = mixtures.SalineWater(network=net)
# Retrieve handles to each species for use below
Cl, Na, H2O = sw.components.values()

# physics
phys = op.physics.GenericPhysics(network=net, phase=sw, geometry=geo)

flow = op.models.physics.hydraulic_conductance.hagen_poiseuille_2D
phys.add_model(propname='throat.hydraulic_conductance',
               pore_viscosity='pore.viscosity',
               throat_viscosity='throat.viscosity',
               model=flow, regen_mode='normal')

current = op.models.physics.ionic_conductance.electroneutrality_2D
phys.add_model(propname='throat.ionic_conductance',
               model=current, regen_mode='normal', ions=[Na.name, Cl.name])

eA_dif = op.models.physics.diffusive_conductance.ordinary_diffusion_2D
phys.add_model(propname='throat.diffusive_conductance.' + Na.name,
               pore_diffusivity='pore.diffusivity.' + Na.name,
               throat_diffusivity='throat.diffusivity.' + Na.name,
               model=eA_dif, regen_mode='normal')

eB_dif = op.models.physics.diffusive_conductance.ordinary_diffusion_2D
phys.add_model(propname='throat.diffusive_conductance.' + Cl.name,
               pore_diffusivity='pore.diffusivity.' + Cl.name,
               throat_diffusivity='throat.diffusivity.' + Cl.name,
               model=eB_dif, regen_mode='normal')

# algorithms
sf = op.algorithms.StokesFlow(network=net, phase=sw)
sf.set_value_BC(pores=net.pores('back'), values=2010)
sf.set_value_BC(pores=net.pores('front'), values=10)
sf.settings['rxn_tolerance'] = 1e-12
sf.run()
sw.update(sf.results())

p = op.algorithms.ChargeConservation(network=net, phase=sw)
p.set_value_BC(pores=net.pores('left'), values=0.02)
p.set_value_BC(pores=net.pores('right'), values=0.01)
p.settings['rxn_tolerance'] = 1e-12
p.settings['charge_conservation'] = 'electroneutrality_2D'

eA = op.algorithms.NernstPlanck(network=net, phase=sw, ion=Na.name)
eA.set_value_BC(pores=net.pores('back'), values=20)
eA.set_value_BC(pores=net.pores('front'), values=10)
eA.settings['rxn_tolerance'] = 1e-12

eB = op.algorithms.NernstPlanck(network=net, phase=sw, ion=Cl.name)
eB.set_value_BC(pores=net.pores('back'), values=20)
eB.set_value_BC(pores=net.pores('front'), values=10)
eB.settings['rxn_tolerance'] = 1e-12

ad_dif_mig_Na = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.' + Na.name,
               model=ad_dif_mig_Na, ion=Na.name,
               s_scheme=scheme)

ad_dif_mig_Cl = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
phys.add_model(propname='throat.ad_dif_mig_conductance.' + Cl.name,
               pore_pressure=sf.settings['quantity'],
               model=ad_dif_mig_Cl, ion=Cl.name,
               s_scheme=scheme)

pnp = op.algorithms.IonicTransport(network=net, phase=sw)
pnp.setup(potential_field=p.name, ions=[eA.name, eB.name])
pnp.settings['i_max_iter'] = 10
pnp.settings['i_tolerance'] = 1e-04

pnp.run()

sw.update(sf.results())
sw.update(p.results())
sw.update(eA.results())
sw.update(eB.results())

# output data
# proj.export_data(phases=[sw], filename='OUT', filetype='xdmf')

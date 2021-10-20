import numpy as np
import openpnm as op


# Workspace and project
ws = op.Workspace()
proj = ws.new_project()

# Network and geometry
np.random.seed(7)
net = op.network.Cubic(shape=[20, 3, 3], spacing=1e-4, project=proj)
geo = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)

# Phase
phase = op.phases.Water(network=net)
phase['pore.consistency'] = 4.2e-2  # Pa.s^n
phase['pore.flow_index'] = 0.52
phase['pore.viscosity_min'] = 0.001
phase['pore.viscosity_max'] = 100

# Physics
phys = op.physics.GenericPhysics(network=net, phase=phase, geometry=geo)

mod1 = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys.add_model(propname='throat.hydraulic_conductance',
               model=mod1, regen_mode='normal')

# Algorithms: Newtonian Stokes flow
sf = op.algorithms.StokesFlow(network=net, phase=phase)
sf.set_value_BC(pores=net.pores('front'), values=1)
sf.set_value_BC(pores=net.pores('back'), values=2)
sf.run()
phase.update(sf.results())
phase['pore.pressure_sf'] = phase['pore.pressure']

mod2 = op.models.physics.hydraulic_conductance.hagen_poiseuille_power_law
phys.add_model(propname='throat.nonNewtonian_hydraulic_conductance',
               model=mod2, regen_mode='normal')

# Algorithms: Non Newtonian Stokes flow
nnsf = op.algorithms.NonNewtonianStokesFlow(network=net, phase=phase)
nnsf.set_value_BC(pores=net.pores('front'), values=1)
nnsf.set_value_BC(pores=net.pores('back'), values=2)
nnsf.run()
phase.update(nnsf.results())

cn = net['throat.conns']
gh_sf = phase['throat.hydraulic_conductance']
gh = phase['throat.nonNewtonian_hydraulic_conductance']
P_sf = phase['pore.pressure_sf']
P = phase['pore.pressure']
Qt_sf = np.abs(gh_sf*np.diff(P_sf[cn], axis=1).squeeze())
Qt = np.abs(gh*np.diff(P[cn], axis=1).squeeze())
Q = Qt/Qt_sf

# Output results to XDMF file format
proj.export_data(phases=[phase], filename='out', filetype='XDMF')

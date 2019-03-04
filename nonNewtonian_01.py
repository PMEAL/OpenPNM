import openpnm as op
import numpy as np

# work space and project
ws = op.Workspace()
proj = ws.new_project()

# network
np.random.seed(7)
net = op.network.Cubic(shape=[20, 1, 1], spacing=1e-4, project=proj)

# geometry
geo = op.geometry.StickAndBall(network=net,
                               pores=net.Ps,
                               throats=net.Ts)

# phase
phase = op.phases.Water(network=net)
phase['pore.consistency'] = 0.0005
phase['pore.flow_index'] = 0.8

# physics
phys = op.physics.GenericPhysics(network=net,
                                 phase=phase,
                                 geometry=geo)

mod1 = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys.add_model(propname='throat.hydraulic_conductance',
               model=mod1, regen_mode='normal')

# algorithms: Newtonian Stokes flow
sf = op.algorithms.StokesFlow(network=net, phase=phase)
sf.set_value_BC(pores=net.pores('front'), values=1)
sf.set_value_BC(pores=net.pores('back'), values=200)
sf.run()
phase.update(sf.results())
phase['pore.pressure_sf'] = phase['pore.pressure']

mod2 = op.models.physics.hydraulic_conductance.hagen_poiseuille_power_law
phys.add_model(propname='throat.nonNewtonian_hydraulic_conductance',
               model=mod2, regen_mode='normal')

# algorithms: Non Newtonian Stokes flow
nnsf = op.algorithms.NonNewtonianStokesFlow(network=net, phase=phase)
nnsf.set_value_BC(pores=net.pores('front'), values=1)
nnsf.set_value_BC(pores=net.pores('back'), values=200)
nnsf.run()
phase.update(nnsf.results())

# output results to a vtk file
# proj.export_data(phases=[phase], filename='OUT', filetype='VTK')

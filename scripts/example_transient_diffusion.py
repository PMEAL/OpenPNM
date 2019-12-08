import openpnm as op
import numpy as np

# work space and project
ws = op.Workspace()
ws.settings["loglevel"] = 30
proj = ws.new_project()

# network
np.random.seed(7)
net = op.network.Cubic(shape=[51, 19, 1], spacing=1e-4, project=proj)

# geometry
geo = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)

# phase
phase = op.phases.Water(network=net)

# physics
phys = op.physics.GenericPhysics(network=net, phase=phase, geometry=geo)
phase['pore.diffusivity'] = 2e-09
phase['throat.diffusivity'] = 2e-09

mod = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance',
               model=mod, regen_mode='normal')

# algorithms: Fickian diffusion
fd = op.algorithms.TransientFickianDiffusion(network=net, phase=phase)
fd.set_value_BC(pores=net.pores('front'), values=0.5)
fd.set_value_BC(pores=net.pores('back'), values=0.1)
fd.setup(t_final=100, t_output=10, t_step=0.5)
fd.run()
phase.update(fd.results())

# output results to a vtk file
phase.update(fd.results())
proj.export_data(phases=[phase], filename='OUT', filetype='xdmf')

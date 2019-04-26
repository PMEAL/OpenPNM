import openpnm as op
import numpy as np

# work space and project
ws = op.Workspace()
ws.settings["loglevel"] = 30
proj = ws.new_project()

# network
np.random.seed(7)
net = op.network.Cubic(shape=[29, 13, 1], spacing=1e-5, project=proj)

# geometry
geo = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)

# phase
phase = op.phases.Water(network=net)

# physics
phys = op.physics.GenericPhysics(network=net, phase=phase, geometry=geo)
phase['pore.diffusivity'] = 2e-09


mod = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance',
               model=mod, regen_mode='normal')

# algorithms: Fickian diffusion
fd = op.algorithms.TransientFickianDiffusion(network=net, phase=phase)
fd.set_value_BC(pores=net.pores('front'), values=0.5)
fd.set_value_BC(pores=net.pores('back'), values=0.2)
fd.set_IC(0.2)
fd.setup(t_scheme='cranknicolson', t_final=100, t_output=5, t_step=1,
         t_tolerance=1e-12)
fd.run()
phase.update(fd.results())

sfd = op.algorithms.FickianDiffusion(network=net, phase=phase)
sfd.set_value_BC(pores=net.pores('front'), values=0.5)
sfd.set_value_BC(pores=net.pores('back'), values=0.2)
sfd.run()

# output results to a vtk file
phase.update(fd.results())
proj.export_data(phases=[phase], filename='OUT', filetype='xdmf')

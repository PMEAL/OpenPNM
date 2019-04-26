import openpnm as op
import numpy as np

# work space and project
ws = op.Workspace()
proj = ws.new_project()

# network
h = 1e-5  # height in m
r = 5e-5  # radius
R = int(round(r/h, 1))  # ratio
nh = 5  # number of pores along the height
nr = int(R*nh)  # number of pores along 1 radius
s = h/(nh-1)  # spacing between pores
cylinder = op.topotools.template_cylinder_annulus(height=nh, outer_radius=nr)
np.random.seed(7)
net = op.network.CubicTemplate(cylinder, spacing=s, project=proj)

Max = net['pore.coords'][:, 2].max()
prs = np.where(net['pore.coords'][:, 2] >= Max)[0]
net['pore.inlet'] = net['pore.all'] * False
net['pore.inlet'][prs] = True

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
c_inlet = 1000  # inlet concentration mol/m3
fd = op.algorithms.TransientFickianDiffusion(network=net, phase=phase)
fd.set_value_BC(pores=net.pores('inlet'), values=c_inlet)
fd.setup(t_final=1, t_output=0.05, t_step=0.05, t_tolerance=1e-12)
fd.run()
phase.update(fd.results())

# output results to a vtk file
phase.update(fd.results())
proj.export_data(phases=[phase], filename='OUT', filetype='xdmf')

import openpnm as op
import numpy as np


# Workspace and project
ws = op.Workspace()
proj = ws.new_project()

# Network
np.random.seed(7)
net = op.network.Cubic(shape=[51, 19, 1], spacing=1e-4, project=proj)

# Geometry
geo = op.geometry.SpheresAndCylinders(network=net, pores=net.Ps, throats=net.Ts)

# Phase
phase = op.phase.Water(network=net)

# Physics
phys = op.physics.GenericPhysics(network=net, phase=phase, geometry=geo)
phase['pore.diffusivity'] = 2e-09
phase['throat.diffusivity'] = 2e-09

mod = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance',
               model=mod, regen_mode='normal')

# Algorithms: Fickian diffusion
fd = op.algorithms.TransientFickianDiffusion(network=net, phase=phase)
fd.set_value_BC(pores=net.pores('front'), values=0.5)
fd.set_value_BC(pores=net.pores('back'), values=0.1)
fd.run(x0=0, tspan=(0, 100), saveat=10)
phase.update(fd.results())

# Output results to a vtk file
phase.update(fd.results())
proj.export_data(phases=[phase], filename='out', filetype='xdmf')

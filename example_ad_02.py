import openpnm as op
import numpy as np

# work space and project
ws = op.Workspace()
proj = ws.new_project()

# network
np.random.seed(7)
n = 10
net = op.network.Cubic(shape=[n, n, 1], spacing=1e-4, project=proj)

# geometry
geo = op.geometry.StickAndBall(network=net,
                               pores=net.Ps,
                               throats=net.Ts)

# phase
phase = op.phases.Water(network=net)

# physics
phys = op.physics.GenericPhysics(network=net,
                                 phase=phase,
                                 geometry=geo)
phase['pore.diffusivity'] = 2e-09
phase['throat.diffusivity'] = 2e-09

mod1 = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys.add_model(propname='throat.hydraulic_conductance',
               model=mod1, regen_mode='normal')

# algorithms: Stokes flow
sf = op.algorithms.StokesFlow(network=net, phase=phase)
sf.set_value_BC(pores=net.pores('left'), values=0)
sf.set_value_BC(pores=net.pores('right'), values=0.1)
sf.run()
phase.update(sf.results())


# algorithms: dispersion
mod2 = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance',
               model=mod2, regen_mode='normal')

mod3 = op.models.physics.ad_dif_conductance.ad_dif
phys.add_model(propname='throat.ad_dif_conductance',
               model=mod3, regen_mode='normal')

dis = op.algorithms.AdvectionDiffusion(network=net, phase=phase)
dis.set_value_BC(pores=net.pores('back'), values=1)
dis.set_value_BC(pores=net.pores('front'), values=0.5)
dis.run()

dis2 = op.algorithms.TransientAdvectionDiffusion(network=net, phase=phase)
dis2.set_value_BC(pores=net.pores('back'), values=1)
dis2.set_value_BC(pores=net.pores('front'), values=0.5)
dis2.setup(t_initial=0, t_final=300, t_step=1, t_output=50, t_tolerance=1e-20,
           t_scheme='implicit')
dis2.run()

# output results to a vtk file
# phase.update(dis.results())
# proj.export_data(phases=[phase], filename='OUT', filetype='VTK')

r"""
Details about the continum and numerical model equations can be found on:
Sadeghi, M. A., Agnaou, M., Barralet, J., & Gostick, J. (2020).

Dispersion modeling in pore networks: A comparison of common pore-scale
models and alternative approaches. Journal of Contaminant Hydrology, 228, 103578.

"""
import openpnm as op
import numpy as np


np.random.seed(0)
ws = op.Workspace()
proj = ws.new_project()

# %% Problem setup
# Create network
net = op.network.Cubic(shape=[33, 33, 1], spacing=9e-4, project=proj)
# Remove the pores at the corners of the network
Ps = (net['pore.back']  * net['pore.right']
    + net['pore.back']  * net['pore.left']
    + net['pore.front'] * net['pore.right']
    + net['pore.front'] * net['pore.left'])
Ts = net['throat.surface']
op.topotools.trim(network=net, pores=net.Ps[Ps], throats=net.Ts[Ts])

# Create geometry
geo = op.geometry.SpheresAndCylinders(network=net, pores=net.Ps, throats=net.Ts)
# Define constant pore/throat diameters
mod = op.models.misc.constant
geo.add_model(propname='pore.diameter', model=mod, value=1.5e-4)
geo.add_model(propname='throat.diameter', model=mod, value=1e-4)
geo.regenerate_models()

# Create phase
water = op.phases.Water(network=net)

# Create physics
phys = op.physics.GenericPhysics(network=net, phase=water, geometry=geo)
# Add model for calculating hydraulic conductance
flow = op.models.physics.hydraulic_conductance.hagen_poiseuille
phys.add_model(propname='throat.hydraulic_conductance',
               pore_viscosity='pore.viscosity',
               throat_viscosity='throat.viscosity',
               model=flow, regen_mode='normal')
# Add model for calculating diffusive conductance
dif = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname='throat.diffusive_conductance', model=dif,
               regen_mode='normal')
# Add model for calculating dispersive conductance (advection-diffusion)
ad_dif = op.models.physics.ad_dif_conductance.ad_dif
phys.add_model(propname='throat.ad_dif_conductance',
               model=ad_dif, regen_mode='normal')
# Add source term
linear = op.models.physics.generic_source_term.linear
phys['pore.A1'] = -1e-15
phys['pore.A2'] = 0.0
phys.add_model(propname='pore.rxn', model=linear, X='pore.concentration',
               A1='pore.A1', A2='pore.A2')
rxn_pores = np.array([200])
net.set_label('rxn', pores=rxn_pores)

# %% Algorithms
# Stokes flow algorithm
sf = op.algorithms.StokesFlow(network=net, phase=water)
sf.set_value_BC(pores=net.pores('back'), values=0.01)
sf.set_value_BC(pores=net.pores('front'), values=0.00)
sf.run()
water.update(sf.results())
# Transient advection diffusion
ad = op.algorithms.TransientAdvectionDiffusion(network=net, phase=water)
ad.set_source(propname='pore.rxn', pores=rxn_pores)
ad.set_value_BC(pores=net.pores('back'), values=100)
ad.set_value_BC(pores=net.pores('front'), values=90)
ad.run(x0=90, tspan=(0, 1000), saveat=200)
water.update(ad.results())

# %% Post processing
# Output results for visualization on Paraview
proj.export_data(phases=[water], filename='out', filetype='xdmf')

# %% Imports
import numpy as np
import openpnm as op
from openpnm.models.physics import source_terms

# %% Initialization: Create Workspace and project objects.
ws = op.Workspace()
ws.settings.default_solver = 'PardisoSpsolve'  # Optionally use ScipySpsolve
ws.settings.loglevel = 50
np.random.seed(9)

# %% Create network, geometry, phase, and physics objects
pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4)
geo = op.geometry.SpheresAndCylinders(network=pn, pores=pn.Ps, throats=pn.Ts)
air = op.phase.Air(network=pn, name='air')
water = op.phase.Water(network=pn, name='h2o')
hg = op.phase.Mercury(network=pn, name='hg')
phys_air = op.physics.Standard(network=pn, phase=air, geometry=geo)
phys_water = op.physics.Standard(network=pn, phase=water, geometry=geo)
phys_hg = op.physics.Standard(network=pn, phase=hg, geometry=geo)

# %% Perform porosimetry simulation
# mip = op.metrics.Porosimetry(network=pn, phase=hg)
# mip.set_inlets(pores=pn.pores(['top', 'bottom']))
# mip.run()
# hg.update(mip.results(Pc=70000))
# mip.plot_intrusion_curve()

# %% Perform Stokes flow simulation
sf = op.algorithms.StokesFlow(network=pn, phase=water)
sf.set_value_BC(pores=pn.pores('left'), values=101325)
sf.set_value_BC(pores=pn.pores('right'), values=0)
sf.run()
water.update(sf.results())
# calculate absolute permeability in x direction
perm = op.metrics.AbsolutePermeability(network=pn)
K = perm.run()
# assert K == 7.51015925e-13

# %% Perform reaction-diffusion simulation
# Add reaction to phys_air
phys_air['pore.n'] = 2
phys_air['pore.A'] = -1e-5
phys_air.add_model(
    propname='pore.2nd_order_rxn',
    model=source_terms.standard_kinetics,
    X='pore.concentration', prefactor='pore.A', exponent='pore.n',
    regen_mode='deferred'
)
# Set up Fickian diffusion simulation
rxn = op.algorithms.FickianDiffusion(network=pn, phase=air)
Ps = pn.find_nearby_pores(pores=50, r=5e-4, flatten=True)
rxn.set_source(propname='pore.2nd_order_rxn', pores=Ps)
rxn.set_value_BC(pores=pn.pores('top'), values=1)
rxn.run()
air.update(rxn.results())

# %% Perform pure diffusion simulation
fd = op.algorithms.FickianDiffusion(network=pn, phase=air)
fd.set_value_BC(pores=pn.pores('left'), values=1)
fd.set_value_BC(pores=pn.pores('right'), values=0)
fd.run()
# calculate formation factor in x direction
FF = op.metrics.FormationFactor(network=pn)
F = FF.run()
# assert F == 20.53387084881872

# %% Output network and the phases to a VTP file for visualization in Paraview
pn.project.export_data(phases=[hg, air, water], filename='output.vtp')

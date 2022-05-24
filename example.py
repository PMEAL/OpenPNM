import openpnm as op
from openpnm.models.physics import source_terms
from openpnm.models import collections
import matplotlib.pyplot as plt
ws = op.Workspace()
ws.clear()


pn = op.network.Cubic(shape=[25, 25, 1])

# Create domain1
Ps = pn.coords[:, 0] < 13
Ts = pn.find_neighbor_throats(pores=Ps, asmask=True)
pn['pore.domain1'] = Ps
pn['throat.domain1'] = Ts

# Create domain2
Ps = pn.coords[:, 0] >= 13
Ts = pn.find_neighbor_throats(pores=Ps, mode='xnor', asmask=True)
pn['pore.domain2'] = Ps
pn['throat.domain2'] = Ts

# Add network/geometry models to both domains
pn.add_model_collection(collections.geometry.cones_and_cylinders, domain='domain1')
pn.add_model_collection(collections.geometry.pyramids_and_cuboids, domain='domain2')

# FIXME: Must regenerate network models, otherwise, phase models will complain
pn.regenerate_models()

# Create phase and add phase/physics models
air = op.phase.Air(network=pn, name="air")
air.models.update(op.models.collections.physics.standard)
air.add_model_collection(collections.physics.standard, domain='domain1')
air.add_model_collection(collections.physics.standard, domain='domain2')
air.regenerate_models()

# Add a nonlinear reaction
air['pore.reaction_sites'] = False
air['pore.reaction_sites'][[310, 212, 113]] = True
air.add_model(propname='pore.reaction',
              model=source_terms.power_law,
              X='pore.concentration',
              domain='reaction_sites',
              A1=-1, A2=2, A3=0, regen_mode='deferred')

# Run Fickian diffusion with reaction
rxn = op.algorithms.FickianDiffusion(network=pn, phase=air)
rxn.set_value_BC(pores=pn.pores('left'), values=1)
rxn.set_source(pores=air.pores('reaction_sites'), propname='pore.reaction')
rxn.run()

# Plot concentration profile
fig, ax = plt.subplots()
ax.pcolormesh(rxn.x.reshape([25, 25]))

import openpnm as op
from openpnm.models.physics import source_terms
from openpnm.models import collections
import matplotlib.pyplot as plt
op.visualization.set_mpl_style()
ws = op.Workspace()
ws.clear()

pn = op.network.Cubic(shape=[25, 25, 1], spacing=1e-4)

# Create domain
Ps = pn.coords[:, 0] < 13e-4
Ts = pn.find_neighbor_throats(pores=Ps, asmask=True)
pn['pore.domain1'] = Ps
pn['throat.domain1'] = Ts

# Create domain2
Ps = pn.coords[:, 0] >= 13e-4
Ts = pn.find_neighbor_throats(pores=Ps, mode='xnor', asmask=True)
pn['pore.domain2'] = Ps
pn['throat.domain2'] = Ts

# Add network/geometry models to both domains
pn.add_model_collection(collections.geometry.cones_and_cylinders,
                        domain='domain1')
pn.add_model_collection(collections.geometry.pyramids_and_cuboids,
                        domain='domain2')

# FIXME: Must regenerate network models, otherwise, phase models will complain
pn.regenerate_models()

# Create phase and add phase/physics models
air = op.phase.Air(network=pn, name="air")
air.add_model_collection(collections.physics.standard)
air.regenerate_models()

# Add a nonlinear reaction
air['pore.reaction_sites'] = False
air['pore.reaction_sites'][[310, 212, 113]] = True
air.add_model(propname='pore.reaction1',
              model=source_terms.power_law,
              X='pore.concentration',
              A1=-1, A2=2, A3=0,
              domain='reaction_sites',
              regen_mode='deferred')
air.add_model(propname='pore.reaction2',
              model=source_terms.power_law,
              X='pore.concentration',
              A1=-1, A2=2, A3=0,
              domain='reaction_sites',
              regen_mode='deferred')

# Run Fickian diffusion with reaction
rxn = op.algorithms.FickianDiffusion(network=pn, phase=air)
rxn.set_value_BC(pores=pn.pores('left'), values=1)
rxn.set_source(pores=air.pores('reaction_sites'), propname='pore.reaction1')
rxn.run()

# Run Fickian diffusion with reaction
rxn2 = op.algorithms.FickianDiffusion(network=pn, phase=air)
rxn2.set_value_BC(pores=pn.pores('left'), values=1)
rxn2.set_source(pores=air.pores('reaction_sites'), propname='pore.reaction2')
rxn2.run()

# Plot concentration profile
fig, ax = plt.subplots()
ax.pcolormesh(rxn.x.reshape([25, 25]))

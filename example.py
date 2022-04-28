import numpy as np
import openpnm as op
import matplotlib.pyplot as plt

pn = op.network.Cubic(shape=[25, 25, 1])
Ps = pn.pores('left')
Ts = pn.find_neighbor_throats(Ps, asmask=True)
pn['throat.left'] = Ts
pn.add_model_collection(op.models.collections.geometry.circles_and_rectangles,
                        domain='all')
pn.regenerate_models()

air = op.phase.Air(network=pn)
air.add_model_collection(op.models.collections.physics.standard, domain='all')
air.regenerate_models()

air['pore.reaction_sites'] = False
air['pore.reaction_sites'][[310, 212, 113]] = True

air.add_model(propname='pore.reaction',
              model=op.models.physics.generic_source_term.power_law,
              X='pore.concentration',
              domain='reaction_sites',
              A1=-1, A2=2, A3=0, regen_mode='deferred')

rxn = op.algorithms.FickianDiffusion(network=pn, phase=air)
rxn.set_value_BC(pores=pn.pores('left'), values=1)
rxn.set_source(pores=air.pores('reaction_sites'), propname='pore.reaction')
rxn.run()
plt.pcolormesh(rxn.x.reshape([25, 25]))

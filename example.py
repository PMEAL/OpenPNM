import numpy as np
import openpnm as op
import matplotlib.pyplot as plt

pn = op.network.Cubic(shape=[25, 25, 1])
pn['pore.domain1'] = pn.coords[:, 0] < 13
pn['throat.domain1'] = pn.find_neighbor_throats(pn.pores('domain1'),
                                                asmask=True)
pn.add_model_collection(op.models.collections.geometry.cones_and_cylinders,
                        domain='domain1')

pn['pore.domain2'] = pn.coords[:, 0] >= 13
pn['throat.domain2'] = pn.find_neighbor_throats(pn.pores('domain2'),
                                                mode='xnor',
                                                asmask=True)
pn.add_model_collection(op.models.collections.geometry.cones_and_cylinders,
                        domain='domain2')

pn.regenerate_models()
# TODO: Need to run twice since there were nans! This needs the dependency
# handler to be fixed up
pn.regenerate_models()

air = op.phase.Air(network=pn)
air.models.update(op.models.collections.physics.standard)
air.add_model_collection(op.models.collections.physics.standard,
                         domain='domain1')
air.add_model_collection(op.models.collections.physics.standard,
                         domain='domain2')
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
# plt.pcolormesh(rxn.x.reshape([25, 25]))

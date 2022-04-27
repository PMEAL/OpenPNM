import numpy as np
import openpnm as op

pn = op.network.Cubic(shape=[4, 4, 1])
Ps = pn.pores('left')
Ts = pn.find_neighbor_throats(Ps, asmask=True)
pn['throat.left'] = Ts
pn.add_model_collection(op.models.collections.geometry.circles_and_rectangles,
                        domain='all')
pn.regenerate_models()

air = op.phase.Air(network=pn)
air.regenerate_models()
air.add_model(propname='throat.diffusive_conductance',
              model=op.models.physics.diffusive_conductance.generic_diffusive)

fd = op.algorithms.FickianDiffusion(network=pn, phase=air)
fd.set_value_BC(pores=pn.pores('left'), values=1)
fd.set_value_BC(pores=pn.pores('right'), values=0)
fd.run()


pn['pore.temperature.air'] = 55
pn['pore.viscosity.air'] = 33
pn['pore.temperature.water'] = 55
pn['pore.viscosity.water'] = 33

a = op.core.Domain()
a.update(pn['*air'])

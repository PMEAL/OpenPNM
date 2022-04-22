import numpy as np
import openpnm as op

pn = op.network.Cubic(shape=[40, 40, 1], name='bob')
pn.models.update(op.models.collections.geometry.circles_and_rectangles)
pn.regenerate_models()

air = op.phase.Air(network=pn)
air.add_model(propname='throat.diffusive_conductance',
              model=op.models.physics.diffusive_conductance.generic_diffusive)
air.run_model('throat.diffusive_conductance')

fd = op.algorithms.FickianDiffusion(network=pn, phase=air)
fd.set_value_BC(pores=pn.pores('left'), values=1)
fd.set_value_BC(pores=pn.pores('right'), values=0)
fd.run()

import matplotlib.pyplot as plt
plt.imshow(fd['pore.concentration'].reshape([40, 40]))

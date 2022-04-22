import numpy as np
import openpnm as op

pn = op.network.Cubic(shape=[4, 4, 1])
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

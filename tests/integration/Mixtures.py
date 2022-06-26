import openpnm as op
from openpnm.models import collections


ws = op.Workspace()
ws.clear()

pn = op.network.Demo(shape=[25, 25, 1], spacing=1e-4)

o2 = op.phase.GasByName(network=pn, species='oxygen', name='O2')
n2 = op.phase.GasByName(network=pn, species='nitrogen', name='N2')

air = op.phase.GasMixture(network=pn, components=[o2, n2], name='air')
air['pore.mole_fraction.O2'] = 0.21
air['pore.mole_fraction.N2'] = 0.79
air.regenerate_models()

air.add_model(propname='throat.diffusive_conductance',
              model=op.models.physics.diffusive_conductance.ordinary_diffusion)

fd = op.algorithms.FickianDiffusion(network=pn, phase=air)
fd.settings['quantity'] = 'pore.concentration.O2'
fd.set_value_BC(pores=pn.pores('left'), values=1.0)
fd.set_value_BC(pores=pn.pores('right'), values=0.0)
fd.run()
print(fd)

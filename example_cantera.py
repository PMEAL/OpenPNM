import openpnm as op
import cantera as ct
from openpnm.utils.misc import tic, toc
ws = op.core.Workspace()

net = op.network.Cubic(shape=[10, 10, 10], spacing=1e-6, name='test_net')
proj = net.project
geom = op.geometry.GenericGeometry(network=net, pores=net.Ps, throats=net.Ts)

air = op.phases.Air(network=net)
phys_air = op.physics.GenericPhysics(network=net, phase=air, geometry=geom)

cantera = op.models.phase.thermal_conductivity.cantera

print('Calculating thermal conductivity with Cantera')
tic()
air.add_model(propname='pore.thermal_conductivity',
              model=cantera,
              cantera_phase_obj=ct.Solution('air.xml'))
air.regenerate_models(propnames='pore.thermal_conductivity')
toc()


print('Calculating thermal conductivity directly')
tic()
water_conductivity = op.models.phase.thermal_conductivity.water
air.add_model(propname='pore.thermal_conductivity', model=water_conductivity)
air.regenerate_models(propnames='pore.thermal_conductivity')
toc()

import numpy as np
import openpnm as op
pn = op.network.Cubic(shape=[10, 10, 10], spacing=0.0001)
Ps = pn.pores('all')
Ts = pn.throats('all')
geo = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts)
hg = op.phases.Mercury(network=pn)
h2o = op.phases.Water(network=pn)
phys1 = op.physics.GenericPhysics(network=pn, phase=hg, geometry=geo)
phys2 = op.physics.GenericPhysics(network=pn, phase=h2o, geometry=geo)
model = op.models.physics.capillary_pressure.washburn
hg.add_model(propname='throat.entry_pressure',
             model=model,
             contact_angle='pore.contact_angle',
             surface_tension='pore.surface_tension')
model = op.models.physics.hydraulic_conductance.hagen_poiseuille
h2o.add_model(propname='throat.hydraulic_conductance',
              model=model,
              throat_viscosity='throat.viscosity',
              pore_diameter='pore.diameter',
              throat_length='throat.length',
              throat_diameter='throat.area')

mip = op.algorithms.Drainage(network=pn)
mip.setup(invading_phase=hg, entry_pressure='throat.entry_pressure')
mip.set_inlets(pn.pores(['left', 'right', 'top', 'bottom', 'front', 'back']))
mip.run(inv_pressures=np.logspace(5, 6.5, num=15))
fig = mip.plot_drainage_curve()
fig.show()

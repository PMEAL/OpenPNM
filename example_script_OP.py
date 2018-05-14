import openpnm as op
ws = op.Workspace()
proj = ws.new_project()
pn = op.network.Cubic(shape=[10, 10, 10], project=proj, spacing=1e-4)
geom = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
geom['pore.volume'][pn.pores('left')] = 0
hg = op.phases.Mercury(network=pn)
phys = op.physics.GenericPhysics(network=pn, phase=hg, geometry=geom)
phys.add_model(propname='throat.entry_pressure',
               model=op.models.physics.capillary_pressure.washburn)
phys.add_model(propname='pore.pc_star',
               model=op.models.misc.from_neighbor_throats,
               throat_prop='throat.entry_pressure',
               mode='min')
phys.add_model(propname='pore.late_filling',
               model=op.models.physics.multiphase.late_filling,
               pressure='pore.pressure',
               Pc_star='pore.pc_star',
               eta=1, Swp_star=0.4)
phys['throat.pc_star'] = phys['throat.entry_pressure']
phys.add_model(propname='throat.late_filling',
               model=op.models.physics.multiphase.late_filling,
               pressure='throat.pressure',
               Pc_star='throat.pc_star',
               eta=1, Swp_star=0.2)

mip = op.algorithms.Porosimetry(project=proj)
mip.set_inlets(pores=pn.pores('left'))
mip.setup(phase=hg)
mip.set_partial_filling(propname='pore.late_filling')
mip.set_partial_filling(propname='throat.late_filling')
mip.run(points=20, stop=1e7)
f = mip.plot_intrusion_curve(fig=f)

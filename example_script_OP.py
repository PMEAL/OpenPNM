import openpnm as op
ws = op.Workspace()
proj = ws.new_project()
pn = op.network.Cubic(shape=[10, 10, 10], project=proj, spacing=1e-4)
geom = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
hg = op.phases.Mercury(network=pn)
phys = op.physics.GenericPhysics(network=pn, phase=hg, geometry=geom)
phys.add_model(propname='throat.capillary_pressure',
               model=op.models.physics.capillary_pressure.washburn)

mip = op.algorithms.SiteAndBondPercolation(project=proj)
mip.set_inlets(pores=pn.pores('left'))
mip.setup(phase=hg,
          access_limited=True,
          mode='bond')
mip.run(points=20, stop=1e8)
mip.plot_percolation_curve()

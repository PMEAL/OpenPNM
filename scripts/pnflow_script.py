import openpnm as op
ws = op.Workspace()
ws.settings['loglevel'] = 30

pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4)
geo = op.geometry._StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
air = op.phases.Air(network=pn)
phys = op.physics.Standard(network=pn, phase=air, geometry=geo)
F = op.algorithms.metrics.FormationFactor(network=pn)
F.run()


op.io.Statoil.add_reservoir_pore(network=pn, pores=pn.pores('left'))
op.io.Statoil.add_reservoir_pore(network=pn, pores=pn.pores('right'))
pn.add_model(propname='throat.total_length',
             model=op.models.geometry.throat_length.ctc)
pn['throat.shape_factor'] = 1.0
pn['pore.shape_factor'] = 1.0
path = "./"
prefix = pn.project.name
op.io.Statoil.export_data(network=pn, path=path, prefix=prefix,
                          shape=[1e-3, 1e-3, 1e-3])
# op.algorithms.metrics.PNFlow.run(pn.project.name, path=path)

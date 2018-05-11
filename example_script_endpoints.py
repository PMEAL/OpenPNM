import openpnm as op
import scipy as sp

ws = op.Workspace()
proj = ws.new_project()
pn = op.network.Cubic(shape=[5, 5, 1], project=proj)
geom = op.geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
geom['pore.diameter'] = sp.rand(pn.Np)
geom.add_model(propname='throat.unit_vector',
               model=op.models.geometry.throat_vector.unit_vector)
geom.add_model(propname='throat.length',
               model=op.models.geometry.throat_length.straight)
geom.add_model(propname='throat.endpoints',
               model=op.models.geometry.throat_length.ends_points,
               throat_vector='throat.unit_vector')
geom.add_model(propname='throat.centroid',
               model=op.models.geometry.throat_vector.centroid)
geom.add_model(propname='throat.conduit_lengths',
               model=op.models.geometry.throat_length.conduit_lengths,
               throat_endpoints='throat.endpoints')
geom.add_model(propname='throat.seed',
               model=op.models.misc.from_neighbor_pores,
               mode='min', pore_prop='pore.diameter')
geom.add_model(propname='throat.diameter',
               model=op.models.misc.scaled,
               factor=0.5, prop='throat.seed')
geom.add_model(propname='throat.area',
               model=op.models.geometry.throat_area.cylinder)
geom.add_model(propname='pore.area',
               model=op.models.geometry.pore_area.sphere)



air = op.phases.Air(network=pn)
phys = op.physics.GenericPhysics(network=pn, phase=air, geometry=geom)
phys.add_model(propname='throat.diffusive_conductance',
               model=op.models.physics.diffusive_conductance.ordinary_diffusion)

alg = op.algorithms.FickianDiffusion(network=pn)
alg.setup(phase=air, quantity='pore.mole_fraction',
          conductance='throat.diffusive_conductance')
alg.set_dirichlet_BC(pores=pn.pores('left'), values=1)
alg.set_dirichlet_BC(pores=pn.pores('right'), values=0)
alg.run()


import matplotlib.pyplot as plt
cmap = plt.cm.jet
norm = plt.Normalize(vmin=sp.amin(pn['pore.diameter']/2),
                     vmax=sp.amax(pn['pore.diameter']/2))
for t in pn.Ts:
    x = geom['throat.endpoints'][t, 0][0], geom['throat.endpoints'][t, 1][0]
    y = geom['throat.endpoints'][t, 0][1], geom['throat.endpoints'][t, 1][1]
    fig = plt.plot(x, y, 'k-o')
for p in pn.Ps:
    xy = pn['pore.coords'][p, :2]
    r = pn['pore.diameter'][p]/2
    circle1 = plt.Circle(xy, r, color=cmap(norm(r)), clip_on=False)
    fig = plt.gcf()
    ax = fig.gca()
    ax.add_artist(circle1)

plt.axis('equal')
plt.axis('off')

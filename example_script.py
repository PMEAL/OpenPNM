import matplotlib.pyplot as plt
import openpnm as op
import scipy as sp
ws = op.Workspace()
proj = ws.new_project()

bp = sp.array([[0.25, 0.25, 0.25], [0.25, 0.75, 0.25],
               [0.75, 0.25, 0.25], [0.75, 0.75, 0.25],
               [0.75, 0.25, 0.75], [0.25, 0.75, 0.75],
               [0.25, 0.25, 0.75], [0.75, 0.75, 0.75]])

s = [1, 1, 1]
#pts = op.topotools.generate_base_points(num_points=200, domain_size=s, reflect=True)
pts = op.topotools.reflect_base_points(base_pts=bp, domain_size=s)

vn = op.network.DelaunayVoronoiDual(points=pts, shape=s)
fig = op.topotools.plot_coordinates(network=vn,
#                                    pores=vn.pores('delaunay'),
                                    color='r')
fig = op.topotools.plot_connections(network=vn,
                                    throats=vn['throat.delaunay'],
                                    color='b',
                                    fig=fig)
fig = op.topotools.plot_connections(network=vn,
                                    throats=vn['throat.voronoi'],
                                    color='g',
                                    fig=fig)
#fig = op.topotools.plot_connections(network=vn,
#                                    throats=vn['throat.interconnect'],
#                                    color='g',
#                                    fig=fig)

#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(xs=pts[:, 0], ys=pts[:, 1], zs=pts[:, 2]).

import openpnm as op
ws = op.Workspace()
proj = ws.new_project()

s = [1]

pts = op.topotools.generate_base_points(num_points=500, domain_size=s, reflect=True)

vn = op.network.DelaunayVoronoiDual(points=pts, shape=s)
fig = op.topotools.plot_connections(network=vn, throats=vn['throat.voronoi'], color='g')
#fig = op.topotools.plot_connections(network=vn, throats=vn['throat.interconnect'], fig=fig, color='r')
#fig = op.topotools.plot_connections(network=vn, throats=vn['throat.delaunay'], fig=fig, color='b')


#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(xs=pts[:, 0], ys=pts[:, 1], zs=pts[:, 2]).


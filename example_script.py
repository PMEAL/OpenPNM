import openpnm as op
ws = op.Workspace()
proj = ws.new_project()

pts = op.topotools.generate_base_points(num_points=1000, domain_size=[1, 1, 1])

dn = op.network.Delaunay(points=pts, shape=[1, 1, 1], trim_domain=True)

#gn = op.network.Gabriel(num_points=100, shape=[1, 1, 1])


vn = op.network.Voronoi(points=pts, shape=[1, 1, 1], trim_domain=True)
op.topotools.plot_connections(network=vn)

dvd = op.network.DelaunayVoronoiDual(points=pts, shape=[1, 1, 1], trim_domain=True)
Ps = dvd.pores('delaunay')
op.topotools.trim(network=dvd, pores=Ps)
op.topotools.plot_connections(network=dvd)

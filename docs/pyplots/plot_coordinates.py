import openpnm as op
pn = op.network.Cubic(shape=[10, 10, 3])
pn.add_boundary_pores()
Ps = pn.pores('internal')
# Create figure showing internal pores
fig = op.topotools.plot_coordinates(network=pn, pores=Ps, c='b')
Ps = pn.pores('*boundary')
# Pass existing fig back into function to plot boundary pores
fig = op.topotools.plot_coordinates(network=pn, pores=Ps, fig=fig, c='r')
fig.show()

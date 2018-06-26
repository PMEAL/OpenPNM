import openpnm as op
pn = op.network.Cubic(shape=[10, 10, 3])
pn.add_boundary_pores()
Ts = pn.throats('*boundary', mode='complement')
# Create figure showing boundary throats
fig = op.topotools.plot_connections(network=pn, throats=Ts)
Ts = pn.throats('*boundary')
# Pass existing fig back into function to plot additional throats
fig = op.topotools.plot_connections(network=pn, throats=Ts, fig=fig, color='r')

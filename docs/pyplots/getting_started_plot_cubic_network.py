import openpnm as op
pn = op.network.Cubic(shape=[10, 10, 10], spacing=0.0001)
fig = op.topotools.plot_connections(network=pn, throats=pn.Ts,
                                    color='b', alpha=0.5)
fig = op.topotools.plot_coordinates(network=pn, pores=pn.Ps, fig=fig,
                                    c='r', s=100)

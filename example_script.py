import openpnm as op

pn = op.network.Cubic(shape=[5, 5, 5])
geom = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

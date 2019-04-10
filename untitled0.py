import openpnm as op
pn = op.network.Cubic(shape=[55, 55, 55], spacing=.0001)
geo = op.geometry.StickAndBall(network=pn)
mip = op.algorithms.metrics.MercuryIntrusion(network=pn)
mip.plot_intrusion_curve()

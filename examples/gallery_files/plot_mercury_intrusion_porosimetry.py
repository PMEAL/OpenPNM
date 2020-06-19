"""
Mercury Intrusion Porosimetry
=============================

This example illustrates in the fewest possible lines how to simulate a mercury intrusion curve on a network

"""

import openpnm as op

pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4)
geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
phase = op.phases.Mercury(network=pn)
phys = op.physics.Standard(network=pn, phase=phase, geometry=geo)
mip = op.algorithms.Porosimetry(network=pn)
mip.setup(phase=phase)
mip.set_inlets(pores=pn.pores('left'))
mip.run()
mip.plot_intrusion_curve()

import openpnm as op
import numpy as np
import matplotlib.pyplot as plt

pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4)
geo = op.geometry.SpheresAndCylinders(network=pn, pores=pn.Ps, throats=pn.Ts)
air = op.phase.Air(network=pn, name='air')
water = op.phase.Water(network=pn, name='h2o')
phys_air = op.physics.Standard(network=pn, phase=air, geometry=geo)
phys_water = op.physics.Standard(network=pn, phase=water, geometry=geo)


ip = op.algorithms.InvasionPercolation(network=pn, phase=water)
ip.set_inlets(pores=pn.pores('left'))
ip.run()


Krel = []
for s in np.linspace(0, pn.Nt, 10):
    inv = ip['throat.invasion_sequence'] < s
    phys_air['throat.hydraulic_conductance'][inv] *= 1e-5
    perm_a = op.algorithms.StokesFlow(network=pn, phase=air)
    perm_a.set_value_BC(pores=pn.pores('top'), values=1)
    perm_a.set_value_BC(pores=pn.pores('bottom'), values=0)
    perm_a.run()
    Krel.append(perm_a.rate(pores=pn.pores('top')))
plt.plot(np.linspace(0, pn.Nt, 10)/pn.Nt, Krel)

# Export to Statoil format.
# Add reservoir pores on each end
op.io.Statoil.add_reservoir_pore(network=pn,
                                 pores=pn.pores('left'),
                                 offset=0.25)
op.io.Statoil.add_reservoir_pore(network=pn,
                                 pores=pn.pores('right'),
                                 offset=0.25)
op.io.Statoil.export_data(network=pn, shape=[10, 10, 10])

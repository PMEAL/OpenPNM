r"""
Example: How to use OpenPNNM to simulating multiphase Fickian diffusion

1D network, the first half of the network is occupied by air and the next
half by water. A partition coefficient of 0.5 is assumed, meaning that the
concentration of the diffusing species in water "at the interface" is half
of that in air.

"""
import openpnm as op
import matplotlib.pyplot as plt
import numpy as np
np.random.seed(10)

# Define network and geometry
net = op.network.Cubic(shape=[10, 1, 1])
geom = op.geometry._StickAndBall(network=net, pores=net.Ps, throats=net.Ts)

# Define constituent phases
air = op.phases.Air(network=net, name="air")
water = op.phases.Water(network=net, name="water")
water["pore.diffusivity"] = air["pore.diffusivity"] * 0.2

# Define MultiPhase object
mphase = op.phases.MultiPhase(network=net, phases=[air, water])
mphase._set_automatic_throat_occupancy()
mphase.set_occupancy(phase=air, pores=[0, 1, 2, 3, 4])
mphase.set_occupancy(phase=water, pores=[5, 6, 7, 8, 9])
const = op.models.misc.constant
mphase.set_binary_partition_coef(phases=[water, air], model=const, value=0.5)

# Define physics object
phys = op.physics.Standard(network=net, phase=mphase, geometry=geom)
mdiff = op.models.physics.diffusive_conductance.multiphase_diffusion
phys.add_model(propname="throat.diffusive_conductance", model=mdiff)

# Define algorithm: Fickian diffusion
fd = op.algorithms.FickianDiffusion(network=net, phase=mphase)
fd.set_value_BC(pores=0, values=1.0)
fd.set_value_BC(pores=9, values=0.1)
fd.run()

# Post-processing
c = fd["pore.concentration"]
plt.figure()
plt.plot(c, "ko:")
plt.xlabel("x (m)")
plt.ylabel("c (mol/m3)")

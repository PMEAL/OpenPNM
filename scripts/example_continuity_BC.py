import openpnm as op
import matplotlib.pyplot as plt
import numpy as np


np.random.seed(10)

# Example 1
net = op.network.Cubic(shape=[10, 1, 1])
phase = op.phases.Air(network=net)
phase["throat.diffusive_conductance"] = 1.0
phase["throat.diffusive_conductance"][4:] = 0.25
phase["throat.diffusive_conductance"][4] = 1e-50
fd = op.algorithms.FickianDiffusion(network=net, phase=phase)
fd.set_value_BC(pores=0, values=1.0)
fd.set_value_BC(pores=9, values=0.0)
fd.set_continuity_BC(ps1=4, ps2=5, K12=2.0)
fd.run()
c = fd["pore.concentration"]
plt.figure()
plt.plot(c, "ko:")
plt.xlabel("x (m)")
plt.ylabel("c (mol/m3)")

# Example 2
net = op.network.Cubic(shape=[10, 10, 1])
xyz = net["pore.coords"]
air = op.phases.Air(network=net)
air["throat.diffusive_conductance"] = 1.0
ts_int = [4, 13, 22, 31, 40, 49, 58, 67, 76, 85]
Ps = xyz[:, 1] < 5
net.set_label("water", pores=Ps)
net.set_label("air", pores=~Ps)
Ts = net.find_neighbor_throats(Ps)
air["throat.diffusive_conductance"][Ts] = 1.0
net.set_label("air_water_int", throats=ts_int)

fd = op.algorithms.FickianDiffusion(network=net, phase=air)
fd.set_value_BC(pores=net.pores(["front", "air"], mode="and"), values=1.0)
fd.set_value_BC(pores=net.pores(["back", "air"], mode="and"), values=0.2)
water_int = net.Ps[(xyz[:, 1] < 5) & (xyz[:, 1] > 4)]
air_int = net.Ps[(xyz[:, 1] > 5) & (xyz[:, 1] < 6)]
net.set_label("water_interface", pores=water_int)
net.set_label("air_interface", pores=air_int)
fd.set_continuity_BC(ps1=air_int, ps2=water_int, K12=2.0)
fd.run()

air.update(fd.results())
c = fd["pore.concentration"]
plt.figure()
plt.imshow(np.rot90(c.reshape(10, 10)))
plt.colorbar()
op.io.XDMF.save(net, phases=air, filename="network")

c1 = c[air_int]
c2 = c[water_int]

A = fd.A.todense()

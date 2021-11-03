import openpnm as op
import numpy as np
import matplotlib.pyplot as plt


# Get Deff w/o including Knudsen effect
spacing = 1.0
net = op.network.Cubic(shape=[10, 10, 10], spacing=spacing)
geom = op.geometry.SpheresAndCylinders(network=net, pores=net.Ps, throats=net.Ts)
air = op.phases.Air(network=net)
phys = op.physics.Standard(network=net, geometry=geom, phase=air)
odiff = op.models.physics.diffusive_conductance.ordinary_diffusion
phys.add_model(propname="throat.diffusive_conductance", model=odiff)
fd = op.algorithms.FickianDiffusion(network=net, phase=air)
fd.set_value_BC(pores=net.pores("left"), values=1.0)
fd.set_value_BC(pores=net.pores("right"), values=0.0)
fd.run()
L = (net.shape * net.spacing)[1]
A = (net.shape * net.spacing)[[0, 2]].prod()
Mdot = fd.rate(pores=net.pores("left")).squeeze()
Deff0 =  Mdot * L / A

# Get Deff w/ including Knudsen effect
mdiff = op.models.physics.diffusive_conductance.mixed_diffusion
phys.add_model(propname="throat.diffusive_conductance", model=mdiff)
spacings = np.linspace(1e-9, 1e-4, 20)
spacings = np.logspace(-9, -3, 25)
Deff = []

for spacing in spacings:
    np.random.seed(10)
    net = op.network.Cubic(shape=[10, 10, 10], spacing=spacing)
    geom = op.geometry.SpheresAndCylinders(network=net, pores=net.Ps, throats=net.Ts)
    air = op.phases.Air(network=net)
    phys = op.physics.Standard(network=net, geometry=geom, phase=air)
    phys.add_model(propname="throat.diffusive_conductance", model=mdiff)
    fd = op.algorithms.FickianDiffusion(network=net, phase=air)
    fd.set_value_BC(pores=net.pores("left"), values=1.0)
    fd.set_value_BC(pores=net.pores("right"), values=0.0)
    fd.run()
    L = (net.shape * net.spacing)[1]
    A = (net.shape * net.spacing)[[0, 2]].prod()
    Mdot = fd.rate(pores=net.pores("left")).squeeze()
    Deff.append(Mdot * L / A)

# Plot ratio of Deff w/ Knudsen to that w/o
Deff = np.array(Deff)
plt.figure()
plt.plot(spacings, Deff/Deff0)
plt.xscale("log")
plt.xlabel("spacing (m)")
plt.ylabel("Deff/Deff0")
